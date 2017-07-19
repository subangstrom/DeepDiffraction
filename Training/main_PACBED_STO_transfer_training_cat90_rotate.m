clear
load('Alex_net.mat')
%demo script for PACBED_rotation CNN transfer training
%note: must use the modified ImageDatastore.m and Trainer.m mode before training
%the modified code is specific for Matlab 2017a
%Weizong Xu, July, 2017
%%
tic
digitDatasetPath_train = 'C:\STEM_pattern_data\STO_PACBED_p227_rotate_13_19mrad';
trainDigitData = imageDatastore(digitDatasetPath_train, 'IncludeSubfolders',true,'LabelSource','foldernames');
toc
%% Make transfer layer
layersTransfer = convnet.Layers(1:end-5);
layersTransfer(end+1) = dropoutLayer(0.5);
layersTransfer(end+1) = convnet.Layers(end-4);
layersTransfer(end+1) = convnet.Layers(end-3);
layersTransfer(end+1) = dropoutLayer(0.5);
layersTransfer(end+1) = fullyConnectedLayer(90,...
              'WeightLearnRateFactor',10,...
	          'BiasLearnRateFactor',20,...
              'name','fc8');

layersTransfer(end+1) = softmaxLayer('name','prob');
layersTransfer(end+1) = classificationLayer('name','classificationLayer');

inputlayer = imageInputLayer([227 227 3],'DataAugmentation',{'randcrop'});
layersTransfer(1)=inputlayer;

optionsTransfer = trainingOptions('sgdm','MaxEpochs',20,...
	'InitialLearnRate',0.0002,...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.5, ...
    'LearnRateDropPeriod', 10,...
    'MiniBatchSize',84,...
    'CheckpointPath',digitDatasetPath_train); 
%% Training
[convnetTransfer, info] = trainNetwork(trainDigitData,layersTransfer,optionsTransfer);
save('try_NN_PACBED_STO_14mrad_Alexnet_cat90_rotate_13-19mrad.mat','convnetTransfer','optionsTransfer','digitDatasetPath_train','info')