classdef Trainer < handle
%The code has been modified for specifit task in the neural network training for Matlab 2017a by Weizong Xu.
%Copy to Matlab installation folder, rename to Trainer.m
%e.g. .../MATLAB_R2017a/toolbox/nnet/cnn/+nnet/+internal/+cnn/Trainer.m

    % Trainer   Class for training a network
    
    %   Copyright 2016 The MathWorks, Inc.
    
    properties(Access = protected)
        Options
        Schedule
        Precision
        Reporter
        ExecutionStrategy
        StopTrainingFlag
    end
    
    methods
        function this = Trainer(opts, precision, reporter, executionSettings)
            % Trainer    Constructor for a network trainer
            %
            % opts - training options (nnet.cnn.TrainingOptionsSGDM)
            % precision - data precision
            this.Options = opts;
            scheduleArguments = iGetArgumentsForScheduleCreation(opts.LearnRateScheduleSettings);
            this.Schedule = nnet.internal.cnn.LearnRateScheduleFactory.create(scheduleArguments{:});
            this.Precision = precision;
            this.Reporter = reporter;
            % Declare execution strategy
            if ismember( executionSettings.executionEnvironment, {'gpu'} )
                this.ExecutionStrategy = nnet.internal.cnn.TrainerGPUStrategy;
            else 
                this.ExecutionStrategy = nnet.internal.cnn.TrainerHostStrategy;
            end
            
            % Print execution environment if in verbose mode
            iPrintExecutionEnvironment(opts, executionSettings);
            
            % Register a listener to detect requests to terminate training
            addlistener( reporter, ...
                'TrainingInterruptEvent', @this.stopTrainingCallback);
        end
        
        function net = train(this, net, data)
            % train   Train a network
            %
            % Inputs
            %    net -- network to train
            %    data -- data encapsulated in a data dispatcher
            % Outputs
            %    net -- trained network
            reporter = this.Reporter;
            schedule = this.Schedule;
            prms = collectSettings(this, net);
            summary = nnet.internal.cnn.util.MiniBatchSummary;
            
            trainingTimer = tic;
            
            reporter.start();
            iteration = 0;
            this.StopTrainingFlag = false;
            [velocity, learnRate] = initializeLearning(this, net);
            this.shuffle( data );
            for epoch = 1:prms.maxEpochs
                data.start();
                while ~data.IsDone && ~this.StopTrainingFlag
                    [X, response] = data.next();
                    % Cast data to appropriate execution environment for
                    % training
                    X = this.ExecutionStrategy.environment(X);
                    X = apply(prms.transforms, X);
                    
                    [gradients, predictions] = this.computeGradients(net, X, response);
                    
                    % Normalize the gradients
                    miniBatchSize = size(X, 4);
                    gradients = this.normalizeGradients(miniBatchSize, gradients);
                    
                    % Reuse the layers outputs to compute loss
                    miniBatchLoss = net.loss( predictions, response );
                    
                    velocity = calculateVelocity( this, ...
                        prms.momentum, velocity, ...
                        prms.l2Regularization, net.LearnableParameters, ...
                        learnRate, gradients);
                    
                    net = net.updateLearnableParameters(velocity);
                    
                    elapsedTime = toc(trainingTimer);
                    
                    iteration = iteration + 1;
                    summary.update(predictions, response, epoch, iteration, elapsedTime, miniBatchLoss, learnRate );
                    reporter.reportIteration( summary );
                end
                learnRate = schedule.update(learnRate, epoch);
                reporter.reportEpoch( epoch, iteration, net );
                info_training=this.Reporter.Reporters{1,3}.Info;
                save('Training_info_temp.mat','info_training')
                % If an interrupt request has been made, break out of the
                % epoch loop
                if this.StopTrainingFlag
                    break;
                end
            end
            reporter.finish();
        end
        
        function net = initializeNetworkNormalizations(this, net, data, precision, executionSettings, verbose)
            
            % Always use 'truncateLast' as we want to process only the data we have.
            savedEndOfEpoch = data.EndOfEpoch;
            data.EndOfEpoch = 'truncateLast';
            
            augmentations = iGetAugmentations(net);
            normalization = iGetNormalization(net);
            
            zerocenter = arrayfun(@(x)isa(x,...
                'nnet.internal.cnn.layer.ZeroCenterImageTransform'), normalization);
            
            if any(zerocenter)
                assert(sum(zerocenter == 1)==1, 'There should only be 1 zero center');
                if verbose
                    iPrintMessage('nnet_cnn:internal:cnn:Trainer:InitializingImageNormalization');
                end
                
                % Weizong modification, allow direct average image input to
                % reduce uncessary waiting time.
                prompt = 'Direct load average image for normalization? Y/N [Y]: ';
                str = input(prompt,'s');
                if isempty(str)
                    str = 'Y';
                end
                
                if strcmp(str,'Y') || strcmp(str,'y') || strcmp(str,'yes') || strcmp(str,'Yes')
                    prompt = 'Load "average_image_data.mat"? Y/N [Y]: ';
                    str = input(prompt,'s');
                    if isempty(str)
                        str = 'Y';
                    end
                    if strcmp(str,'Y') || strcmp(str,'y') || strcmp(str,'yes') || strcmp(str,'Yes')
                        load('average_image_data.mat','avgI');
                        avgI=double(avgI);
                    else
                        chk_err=1;
                        while chk_err==1
                            prompt = 'Please input file name (include .mat) ';
                            filename = input(prompt,'s');
                            try
                                load(filename,'avgI');
                                avgI=double(avgI);
                                chk_err=0;
                            catch ME
                                disp(ME)
                                chk_err=1;
                            end
                        end
                    end
                else
                    avgI = this.ExecutionStrategy.computeAverageImage(data, augmentations, executionSettings);
                end
                %  avgI = this.ExecutionStrategy.computeAverageImage(data, augmentations, executionSettings);
                option_save=this.Options;
                save(['average_image_data_',datestr(now,'yyyymmdd_HH_MM_SS'),'.mat'],'avgI','augmentations','option_save');
                net.Layers{1}.AverageImage = precision.cast(avgI);
            end
            
            data.EndOfEpoch = savedEndOfEpoch;
        end
    end
    
    methods(Access = protected)
        function stopTrainingCallback(this, ~, ~)
            % stopTraining  Callback triggered by interrupt events that
            % want to request training to stop
            this.StopTrainingFlag = true;
        end
            
        function settings = collectSettings(this, net)
            % collectSettings  Collect together fixed settings from the
            % Trainer and the data and put in the correct form.
            settings.maxEpochs = this.Options.MaxEpochs;
            settings.lossFunctionType = iGetLossFunctionType(net);
            settings.shuffleOption = this.Options.Shuffle;
            
            settings.momentum = this.Precision.cast( this.Options.Momentum );
            settings.l2Regularization = this.Precision.cast( this.Options.L2Regularization );
            try
                settings.transforms = [iGetAugmentations(net) iGetNormalization(net)];
            catch
                settings.transforms = [iGetAugmentations(net)' iGetNormalization(net)];
            end
        end
        
        function [velocity, learnRate] = initializeLearning(this, net)
            % initializeLearning  Set the learning parameters to their
            % starting values.
            velocity = iInitializeVelocity(net, this.Precision);
            learnRate = this.Precision.cast( this.Options.InitialLearnRate );
        end
        
        function [gradients, predictions] = computeGradients(~, net, X, Y)
            % computeGradients   Compute the gradients of the network. This
            % function returns also the network output so that we will not
            % need to perform the forward propagation step again.
            [gradients, predictions] = net.computeGradientsForTraining(X, Y);
        end
        
        function gradients = normalizeGradients(~, miniBatchSize, gradients)
            % normalizeGradients   Normalize the gradients by dividing by
            % number of examples in the mini batch
            gradients = cellfun(@(x)x/miniBatchSize, gradients, 'UniformOutput', false);
        end
        
        function newVelocity = calculateVelocity(~, momentum, oldVelocity, globalL2Regularization, learnableParametersArray, globalLearnRate, gradients)
            numLearnableParameters = numel(learnableParametersArray);            
            newVelocity = cell(numLearnableParameters, 1);
            for i = 1:numLearnableParameters
                param = learnableParametersArray(i);
                newVelocity{i} = iVelocity(momentum,oldVelocity{i}, ...
                    globalL2Regularization, globalLearnRate, gradients{i}, ...
                    param.L2Factor, param.LearnRateFactor, param.Value);
            end
        end
    end
    
    methods(Access = private)
        function shuffle(this,data)
            % shuffle   Shuffle the data at start of training as per
            % training options
            if isequal(this.Options.Shuffle, 'once')
                data.shuffle();
            end
        end
    end
end

function t = iGetLossFunctionType(net)
if isempty(net.Layers)
    t = 'nnet.internal.cnn.layer.NullLayer';
else
    t = class(net.Layers{end});
end
end

function n = iGetNormalization(net)
if isempty(net.Layers)
    n = nnet.internal.cnn.layer.ImageTransform.empty;
else
    n = net.Layers{1}.Transforms;
end
end

function a = iGetAugmentations(net)
if isempty(net.Layers)
    a = nnet.internal.cnn.layer.ImageTransform.empty;
else
    a = net.Layers{1}.TrainTransforms;
end
end

function scheduleArguments = iGetArgumentsForScheduleCreation(learnRateScheduleSettings)
scheduleArguments = struct2cell(learnRateScheduleSettings);
end

function velocity = iInitializeVelocity(net, precision)
velocity = num2cell( precision.cast(zeros(numel(net.LearnableParameters),1)) );
end

function vNew = iVelocity(m,vOld,gL2,gLR,grad,lL2,lLR,W)
% [1]   A. Krizhevsky, I. Sutskever, and G. E. Hinton, "ImageNet
%       Classification with Deep Convolutional Neural Networks", in
%       Advances in Neural Information Processing Systems 25, 2012.

% m = momentum
% gL2 = globalL2Regularization
% gLR = globalLearnRate
% g = gradients, i.e., deriv of loss wrt weights
% lL2 = l2Factors, i.e., learn rate for a particular factor
% lLR = learnRateFactors
% W = learnableParameters

% learn rate for this parameters
alpha = gLR*lLR;
% L2 regularization for this parameters
lambda = gL2*lL2;

% Velocity formula as per [1]
vNew = m*vOld - lambda*alpha*W - alpha*grad;
end

function iPrintMessage(messageID, varargin)
string = getString(message(messageID, varargin{:}));
fprintf( '%s\n', string );
end

function iPrintExecutionEnvironment(opts, executionSettings)
% Print execution environment if in 'auto' mode
if opts.Verbose 
    if ismember(opts.ExecutionEnvironment, {'auto'})
        if ismember(executionSettings.executionEnvironment, {'cpu'})
            iPrintMessage( ...
                'nnet_cnn:internal:cnn:Trainer:TrainingInSerialOnCPU');
        elseif ismember(executionSettings.executionEnvironment, {'gpu'})
            iPrintMessage( ...
                'nnet_cnn:internal:cnn:Trainer:TrainingInSerialOnGPU');
        end
    elseif ismember(opts.ExecutionEnvironment, {'parallel'})
        if ismember(executionSettings.executionEnvironment, {'cpu'})
            iPrintMessage( ...
                'nnet_cnn:internal:cnn:Trainer:TrainingInParallelOnCPUs');
        elseif ismember(executionSettings.executionEnvironment, {'gpu'})
            iPrintMessage( ...
                'nnet_cnn:internal:cnn:Trainer:TrainingInParallelOnGPUs');
        end
    end
end
end
