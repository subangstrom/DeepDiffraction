%Hybrid LSF and CNN method
%Weizong Xu, wxu4@ncsu.edu, July, 2017
clear
%Matlab 2017a and later version is recommended.
%Please add subfolder to path.

%% LSF Parameter Setup (to get initial range for LSF)
LSF_thickness_search=6; %+- nm thickness range for the search
search_range=[1,120]; %search range for thickness
rot_angle_list=-1:1:1; %rotation angle offset during search
img_size_list=0.680:0.01:0.720; %size search range for opt size (CNN aligned at around 320-330 size)
int_comp_list=0.9:0.075:1.275; %intensity scale range for search
offsetx_list=-2:1:2; %shift x offset during search
offsety_list=-2:1:2;%shift y offset during search
crystal_symmetry=4; %support 4:four-fold, 2:two-fold

%% CNN - load NN model for PACBED and database
% net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_19D1mrad.mat');%for 19.1mrad SrTiO3
% net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_19D1mrad.mat');%for 19.1mrad SrTiO3
net_str.convnet_thickness=CNN_net_load('PACBED_ConvNet_thickness_13D6mrad.mat');%for 13.6mrad SrTiO3
net_str.convnet_tilt=CNN_net_load('PACBED_ConvNet_tilt_13D6mrad.mat');%for 13.6mrad SrTiO3
net_str.convnet_shift=CNN_net_load('PACBED_ConvNet_shift.mat');
net_str.convnet_size=CNN_net_load('PACBED_ConvNet_size.mat');
net_str.convnet_rotate=CNN_net_load('PACBED_ConvNet_rotate.mat');

%% load database
load('Database_14mrad_STO_demo.mat','database','name_list')
% load('Database_19mrad_STO_demo.mat','database','name_list')

%% CNN - Image alignment and thickness/tilt measurement from CNN
option.mode='process_batch';
option.shift_enhance=1; %1-on, others off
option.size_stable=0; %if all PACBEDs are in same size (2D map mode)
option.rotation_stable=0; % output PACBEDs in same rotation angle (2D map mode)
option.shift_avg=0; %average shift using 2D neighbour info
option.crop_size_opt=330;
option.output_rot_angle=0; %reset output PACBED angle to value set, by default nan, i.e. not reset angle 
option.iter_num_max=6; %4-6 good enough
[ img_align_cell, para_in, option ] = align_PACBED( database, net_str, option );
option.thickness_comp_ratio=0.70;
option.angle_list=0:30:60;
option.size_list=-5:5:5;
option.comp_size_thickness=-5; %(-10%~0) is better for CNN thickness(%)
option.comp_size_tilt=0; %(0-5%) is better for CNN tlit(%)
option.symmetry=4; %STO [001] 4:four-fold; STO [011] 2:two-fold
[ thickness_determine_CNN, tilt_determine_CNN, data_out_all_CNN, img_align_cell_tilt] = database_identify_enhance_batch( database, net_str, para_in, option );

%% CNN - convert CNN output for LSF code
size_input=size(database);
data_filename=name_list;
thickness_determine=thickness_determine_CNN(:,1);
tilt_determine_H=tilt_determine_CNN(:,1);
tilt_determine_G=tilt_determine_CNN(:,2);
tilt_determine_P=tilt_determine_CNN(:,5);
chk_mark=tilt_determine_CNN(:,6);
tilt_determine=cell(length(tilt_determine_H),2);
for i=1:length(tilt_determine_H) %format convert
    tilt_determine{i,1}(1)=min(abs([tilt_determine_H(i,1),tilt_determine_G(i,1)]));
    tilt_determine{i,1}(2)=max(abs([tilt_determine_H(i,1),tilt_determine_G(i,1)]));
    tilt_determine{i,2}=tilt_determine_P(i,1)/100;
end

%% LSF - get tilt label number and its neighbour list
load('PACBED_STO_13D6mrad_data_1nm_resolution_mrad_tilt_38cat.mat','PACBED_data','name_list','tilt_angle','tilt_nei_list_opt','tilt_nei_list')
for i=1:size(tilt_determine,1)
    for j=1:length(tilt_angle)
        if tilt_determine{i,1}(1)==tilt_angle{j}(1) && tilt_determine{i,1}(2)==tilt_angle{j}(2)
            tilt_determine{i,3}=j;
            tilt_determine{i,4}=tilt_nei_list_opt{j}(:,2)'; %optimal search range
%             tilt_determine{i,4}=tilt_nei_list{j}(:,2)'; %full range mistilt search
        end
    end
end
%% LSF excess background remove (high pass filter)
[ img_align_cell_tilt, background ] = image_series_BG_rm( img_align_cell_tilt, 0 ); %0 simple background substract on corner intensity
Img_series_Looper(img_align_cell_tilt,data_filename)
% Img_series_Looper(background,data_filename)

%% LSF - main
data_out=cell(length(data_filename),5);time_count=[];
for num_select=1:length(data_filename)
%     tic
    filename=data_filename{num_select};
    disp(['Load processed image from ',filename])
    A=img_align_cell_tilt{num_select,1}(:,:,1);
    thickness_est=thickness_determine(num_select,1);
    comp_int=uint8(floor(thickness_est*0.4));%substract background according to the determined thickness
    img_in=A-comp_int;

    s_list.name_list=name_list;
    s_list.filename=filename;
    s_list.chk_mark=chk_mark(num_select); %chkmark 0:CNN pre-determined, 1:CNN uncertain, 2: full range (4-fold symmetry) -1: full range (2-fold symmetry) 
    if crystal_symmetry~=4 %|| crystal_symmetry~=2
        s_list.chk_mark=-1; %if symmettry is not the current support version of trained CNN
    end
    %first try
    t_select_list=max(thickness_est-LSF_thickness_search,min(search_range)):2:min(thickness_est+LSF_thickness_search,max(search_range));
    tilt_select_list=tilt_determine{num_select,4};%1:size(PACBED_data,1); get closer neighbour list
    s_list.rot_angle_list=rot_angle_list;
    s_list.img_size_list=img_size_list;
    s_list.int_comp_list=int_comp_list;
    s_list.tilt_select_list=tilt_select_list;
    s_list.t_select_list=t_select_list;
    s_list.offsetx_list=offsetx_list;
	s_list.offsety_list=offsety_list;
    threshold=1.4e7;
    %1st try
    tic
    [ data_out_try ] = func_lsf_PACBED( img_in, PACBED_data, s_list, threshold );

    %get closest range
    threshold=data_out_try{1,1}(1,1);
    s_list.rot_angle_list=data_out_try{1,1}(1,2)-1:1:data_out_try{1,1}(1,2)+1;
    s_list.img_size_list=data_out_try{1,1}(1,3)-0.01:0.0025:data_out_try{1,1}(1,3)+0.01;
    s_list.int_comp_list=data_out_try{1,1}(1,4)-0.1:0.025:data_out_try{1,1}(1,4)+0.1;
    s_list.t_select_list=max(data_out_try{1,1}(1,5)-2,min(search_range)):1:min(data_out_try{1,1}(1,5)+2,max(search_range));
    s_list.offsetx_list=data_out_try{1,1}(1,6)-1:1:data_out_try{1,1}(1,6)+1;
    s_list.offsety_list=data_out_try{1,1}(1,7)-1:1:data_out_try{1,1}(1,7)+1;

    [ data_out_single ] = func_lsf_PACBED( img_in, PACBED_data, s_list, threshold );
    toc
    time_count(num_select)=toc;

    AB_list_sort=data_out_single{1,1};
    t_select=AB_list_sort(1,5);
    tilt_select=AB_list_sort(1,8);
    disp(['It is ', name_list{tilt_select,t_select}])
    disp('---------------------------')
    disp(' ')
    data_out(num_select,:)=data_out_single;
end
thickness_determine_LSF=cell2mat(data_out(:,4));
% save('least_square_fit_1D_14mrad_hybrid.mat','data_out')

%% tilt analysis
tilt_determine_LSF=data_out(:,3);
for i=1:length(tilt_determine_LSF)
    a=tilt_determine_LSF{i};
    num1=regexp(a, '_');
    a=a(num1(1)+1:num1(2)-1);
    num=regexp(a, 'G');
    str1=a(2:num-1);
    str1=strrep(str1,'D','.');
    str2=a(num+1:end);
    str2=strrep(str2,'D','.');
    tilt_determine_LSF{i,1}=[str2double(str1),str2double(str2)];
end
[ tilt_HG, tilt_r, tilt_azimuth, img_align_cell_tilt2 ] = tilt_convert_full_coor( img_align_cell, tilt_determine_LSF, crystal_symmetry);

%% Look at the best match
chk_img_database={};chk_name_list={};
for num_select=1:size(data_out,1)
    LSF_opt=data_out{num_select,1};
    img_in=img_align_cell_tilt2{num_select,1}(:,:,1);
    [ A2, B_crop, AB_diff ] = func_lsf_PACBED_evaluate( img_in, PACBED_data, LSF_opt, crystal_symmetry );
    chk_img_database{1,num_select}=A2;
    chk_name_list{1,num_select}=data_out{num_select,2};
    chk_img_database{2,num_select}=B_crop;
    chk_name_list{2,num_select}=data_out{num_select,3};
    chk_img_database{3,num_select}=AB_diff;
    chk_name_list{3,num_select}='Difference';
end
Img_series_Looper(chk_img_database,chk_name_list)