clear
load('Alex_net.mat')
%demo script for PACBED_tilt(13.6mrad) CNN transfer training
%note: must use the modified ImageDatastore.m and Trainer.m mode before
%training, see Tips for training.pdf
%the modified code is specific for Matlab 2017a
%Weizong Xu, July, 2017
%% Search training data and generate data list
tic
digitDatasetPath_train = 'D:\STEM_pattern_data\STO_PACBED_p227_AlexNet_t_cat38_4mrad';
trainDigitData = imageDatastore(digitDatasetPath_train, 'IncludeSubfolders',true,'LabelSource','foldernames');
toc
%% Make transfer layer
layersTransfer = convnet.Layers(1:end-5);%retain these layers from AlexNet
layersTransfer(end+1) = dropoutLayer(0.5); %insert dropout layer after fully connected layer to reduce overfitting during training
layersTransfer(end+1) = convnet.Layers(end-4);% add fully connected layer from AlexNet 
layersTransfer(end+1) = convnet.Layers(end-3);% add fully connected layer from AlexNet 
layersTransfer(end+1) = dropoutLayer(0.5);%insert dropout layer after fully connected layer to reduce overfitting during training
layersTransfer(end+1) = fullyConnectedLayer(38,... % add fully connected layer based on the total number of categories.
              'WeightLearnRateFactor',10,... %10x faster than other base layers
	          'BiasLearnRateFactor',20,...   %20x faster than other base layers
              'name','fc8');

layersTransfer(end+1) = softmaxLayer('name','prob'); %add softmax for classification
layersTransfer(end+1) = classificationLayer('name','classificationLayer');

inputlayer = imageInputLayer([227 227 3],'DataAugmentation',{'randfliplr','randcrop'});%data augumentation for shift and crop
layersTransfer(1)=inputlayer;

optionsTransfer = trainingOptions('sgdm','MaxEpochs',10,... %max 10 epochs
	'InitialLearnRate',0.00005,...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 10,...
    'MiniBatchSize',84,... %opt at 84 for fastest speed
    'CheckpointPath',digitDatasetPath_train); %output data after every 1 epoch
%% Training
[convnetTransfer, info] = trainNetwork(trainDigitData,layersTransfer,optionsTransfer);
save('try_NN_PACBED_STO_14mrad_Alexnet_cat38_tilt_4mrad.mat','convnetTransfer','optionsTransfer','digitDatasetPath_train','info')