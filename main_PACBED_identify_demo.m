%CNN measurement for PACBED (parallel edition)
%Weizong Xu, wxu4@ncsu.edu, July, 2017
clear
%Matlab 2017a and later version is recommended.
%Please add subfolder to path before run. Run section by section is recommended.
%Please contact Prof. James LeBeau (jmlebeau@ncsu.edu) for large neural network files.
%% load NN model for PACBED
net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_13D6mrad.mat');%for 13.6mrad SrTiO3
net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_13D6mrad.mat');%for 13.6mrad SrTiO3
% net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_19D1mrad.mat');%for 19.1mrad SrTiO3
% net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_19D1mrad.mat');%for 19.1mrad SrTiO3
net_str.convnet_shift=CNN_net_load('PACBED_ConvNet_shift.mat');
net_str.convnet_size=CNN_net_load('PACBED_ConvNet_size.mat');
net_str.convnet_rotate=CNN_net_load('PACBED_ConvNet_rotate.mat');

%% load data
load('Database_14mrad_STO_demo.mat','database','name_list')
% load('Database_19mrad_STO_demo.mat','database','name_list')
% [database, name_list] = file_load; %or load image files manually
Img_series_Looper(database',name_list')

%% PACBED alignment (speed-up version)
option.mode='process_batch';
option.cal_gpu=1; %1-gpu others cpu
option.shift_enhance=1; %1-on, others off
option.size_stable=0; %if all PACBEDs are in same size (2D map mode)
option.rotation_stable=0; % output PACBEDs in same rotation angle (2D map mode)
option.shift_avg=0; %average shift using 2D neighbour info
option.crop_size_opt=330; %optimal size 330
option.iter_num_max=6; %4-6 good enough
option.output_rot_angle=0; %reset output PACBED angle to value set, by default nan, i.e. not reset angle 
[ img_out, para_in, option ] = align_PACBED( database, net_str, option );
Img_series_Looper(img_out,name_list);
%% thickness/tilt measurment
option.thickness_comp_ratio=0.7; %rm dI=ratio*thickness - 0.7 optimized for 13.6mrad, 0.3 for 19.1mrad
option.angle_list=0:30:60; %angle pooling for more robust result, only 0 for faster process
option.size_list=-5:5:5; %size pooling for more robust result, only 0 for faster process
option.thickness_stable=0; %+- value around median is accepted, 0-disable this funtion
[ thickness_determine, tilt_determine] = database_identify_enhance_batch( database, net_str, para_in, option );
%% Display thickness plot (1D) if make sense from a line scan data
figure;plot(thickness_determine(:,1));title('Thickness');
figure;plot(tilt_determine(:,1));title('H value');
figure;plot(tilt_determine(:,2));title('G value');
figure;plot(tilt_determine(:,3));title('Tilt amplitude(mrad)');
figure;plot(tilt_determine(:,4));title('azimuth angle (degree)')