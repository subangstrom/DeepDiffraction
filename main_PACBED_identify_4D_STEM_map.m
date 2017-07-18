%CNN measurement for PACBED (parallel edition for 4D STEM data)
%Weizong Xu, wxu4@ncsu.edu, July, 2017
clear
%Matlab 2017a and later version is recommended.
%Please add subfolder to path.

%% load NN model for PACBED
% net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_19D1mrad.mat');%for 19.1mrad SrTiO3
% net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_19D1mrad.mat');%for 19.1mrad SrTiO3
net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_13D6mrad.mat');%for 13.6mrad SrTiO3
net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_13D6mrad.mat');%for 13.6mrad SrTiO3
net_str.convnet_shift=CNN_net_load('PACBED_ConvNet_shift.mat');
net_str.convnet_size=CNN_net_load('PACBED_ConvNet_size.mat');
net_str.convnet_rotate=CNN_net_load('PACBED_ConvNet_rotate.mat');

%% load database
load('4D-STEM_area_demo.mat','database','name_list')
Img_series_Looper(database',name_list')

%% PACBED alignment (speed-up version)
option.mode='process_batch';
option.shift_enhance=1; %1-on, others off
option.size_stable=1; %if all PACBEDs are in same size (map mode)
option.rotation_stable=1; % output PACBEDs in same rotation angle (2D map mode)
option.shift_avg=3; %average shift using 2D neighbour info
option.crop_size_opt=330; %optimized size is around 320-330
[ img_align_cell, para_in, option ] = align_PACBED( database, net_str, option );
Img_series_Looper(img_align_cell,name_list);

%% thickness/tilt measurment
option.thickness_comp_ratio=0.70; %rm dI=ratio*thickness
option.angle_list=0:30:90; %angle pooling for more robust result, only 0 for faster process
option.size_list=-5:5:5; %size pooling for more robust result, only 0 for faster process
option.thickness_stable=0; %+- value around median is accepted, 0-disable this funtion
option.comp_size_thickness=-5; %size compensation, %(-10%~0) is better for CNN thickness
option.comp_size_tilt=0; %size compensation, %(0~5) is better for CNN tilt
[ thickness_determine, tilt_determine] = database_identify_enhance_batch( database, net_str, para_in, option );

%% Display thickness map
figure;imagesc(thickness_determine(:,:,1));axis image;title('thickness (nm)');colorbar;
figure;imagesc(tilt_determine(:,:,1));axis image;title('H value');colorbar;
figure;imagesc(tilt_determine(:,:,2));axis image;title('G value');colorbar;
figure;imagesc(tilt_determine(:,:,3));axis image;title('Tilt amplitude(mrad)');colorbar;
figure;imagesc(tilt_determine(:,:,4));axis image;title('azimuth angle (degree)');colorbar;