clear
load('Alex_net.mat')
%demo script for PACBED_thickness(19.1mrad) CNN transfer training
%note: must use the modified ImageDatastore.m and Trainer.m mode before training
%the modified code is specific for Matlab 2017a
%Weizong Xu, July, 2017
%%
tic
digitDatasetPath_train = 'E:\STEM_pattern_data\STO_PACBED_p227_19mrad_thickness';
trainDigitData = imageDatastore(digitDatasetPath_train, 'IncludeSubfolders',true,'LabelSource','foldernames');
toc
%% Make transfer layer
layersTransfer = convnet.Layers(1:end-5);
layersTransfer(end+1) = dropoutLayer(0.5);
layersTransfer(end+1) = convnet.Layers(end-4);
layersTransfer(end+1) = convnet.Layers(end-3);
layersTransfer(end+1) = dropoutLayer(0.5);
layersTransfer(end+1) = fullyConnectedLayer(80,...
              'WeightLearnRateFactor',10,...
	          'BiasLearnRateFactor',20,...
              'name','fc8');

layersTransfer(end+1) = softmaxLayer('name','prob');
layersTransfer(end+1) = classificationLayer('name','classificationLayer');

inputlayer = imageInputLayer([227 227 3],'DataAugmentation',{'randfliplr','randcrop'});
layersTransfer(1)=inputlayer;

optionsTransfer = trainingOptions('sgdm','MaxEpochs',10,...
	'InitialLearnRate',0.000005,...
    'MiniBatchSize',84,...
    'CheckpointPath',digitDatasetPath_train); 
%% Training
[convnetTransfer, info] = trainNetwork(trainDigitData,layersTransfer,optionsTransfer);
save('try_NN_PACBED_STO_19D1_mrad_Alexnet_cat80_thickness.mat','convnetTransfer','optionsTransfer','digitDatasetPath_train','info')