clear; %close all;clc
%demo script for PACBED_shift CNN data prep
%Weizong Xu, July, 2017
%Please contact Prof. James LeBeau (jmlebeau@ncsu.edu) for simulation files.
%% get more images (generate signal images)
load('PACBED_STO_14mrad_data_1nm_resolution_0-4mrad_tilt_corr3.mat')
for i=1:size(PACBED_data,1) %reduce image size to save memory
    for j=1:size(PACBED_data,2)
        if size(PACBED_data{i,j},1)==1032 %938 - 14.5mrad
            PACBED_data{i,j}=PACBED_data{i,j}(128:end-127,128:end-127,1);
        else
            PACBED_data{i,j}=PACBED_data{i,j}(81:end-80,81:end-80,1);
        end
    end
end
tmp_a=PACBED_data;
tmp_b=name_list;
tmp_c=tilt_list';
load('PACBED_STO_19mrad_data_1nm_resolution_0-4mrad_tilt_corr3.mat')
PACBED_data=[PACBED_data; tmp_a];
name_list=[name_list; tmp_b];
tilt_list=[tilt_list'; tmp_c];
clear tmp_a tmp_b tmp_c
%%
img_num=1; %number of noisy image in the same condition
fpwd='F:\STEM_pattern_data\STO_PACBED_p227_shift_14_19mrad';%Folder of data for training
num_center=15;%number of random images per condition
index_list=[0];%noise level from 0-25 (0-no noise, 25-noisest)
shift_list=-30:2:30;%set shift categories
for i=1:length(shift_list)
    if ~exist([fpwd,'\','STO_shift_',num2str(shift_list(i))],'dir')
        disp([fpwd,'\','STO_shift_',num2str(shift_list(i)),' is created.'])
    end
end
%%
tic
for i_shift=1:length(shift_list) %parfor compatible, but it is recommended to run multiple instances
    shift_x=shift_list(i_shift);
    image_sum=zeros(227,227,3);
    num=uint32(0);
    for i_size=240:20:430 %set size range, value is the crop image pixel size
        for i_list=6:2:length(list_total) %set thickness range, only 2 nm resolution is good enough
            for i_tilt=1:size(PACBED_data,1)
                PACBED_read=PACBED_data{i_tilt,i_list};
                option=[];
                option.crop_size=i_size;
                rng('shuffle')
                option.rng_set=0; % run shuffle gloably to speed up
                option.noise_type='gaussian';%'poisson' is also supported;
                option.rot_Angle=0;
                option.image_output_size=227;
                option.image_scale=255;
                option.noise_level=0;
                option.signal_level=10;
                option.chk_print=0;
                option.Gaus_blur_size=1;
                option.Avg_size=1;
                option.image_shift=[0,0];
                image_size=size(PACBED_read,1);
                
                i_inten=1.0; %set 1, i.e. no intensity scalling, training programm will do such job
                option.max_intensity_ratio=i_inten;
                disp(['Shift',num2str(shift_x),'Crop size:',num2str(i_size),' ',num2str(i_list),'nm ',tilt_list{i_tilt},'_TiltCat_',num2str(i_inten), ' #', num2str(num)])
                for i_index=1:length(index_list)
                    option.noise_index=index_list(i_index);
                    for i_Center=1:num_center
                        img_center=[round(image_size/2), round(image_size/2)];
                        rand_x=image_size;rand_y=image_size;
                        while shift_x*shift_x+rand_y*rand_y>2500 %set random shift range
                            rand_y=round((rand-0.5)*80);
                            option.signal_center(1)=shift_x+img_center(1);
                            option.signal_center(2)=rand_y+img_center(2);
                        end
                        option.distort_signal=1;
                        option.distort_x=(rand()-0.5)*0.02;%set random distorion via affine transformation
                        option.distort_y=(rand()-0.5)*0.02;
                        if i_Center==1 %be sure one training data is in the center
                            option.signal_center(2)=img_center(2);
                            option.distort_x=0;
                            option.distort_y=0;
                        end
                        i_angle=round(rand()*360);% set rotation angle randomly
                        PACBED_read_rot=imrotate(PACBED_read,i_angle,'bicubic','crop');
                        [ image_out, image_raw, ~ ] = PACBED_add_noise_gray( PACBED_read_rot, img_num, option );%add noise, distortion to the image
                        
                        if i_index>0
                            img_num_eff=img_num;
                        else
                            img_num_eff=1;
                        end
                        %output
                        for i_img=1:img_num_eff
                            num=num+1;
                            image_sum=image_sum+double(image_out{i_img});
                            filename=[fpwd,'/','STO_shift_',num2str(shift_list(i_shift)),'/Img_',num2str(num),'.jpg'];
                            imwrite(uint8(image_out{i_img}),filename, 'Quality',94) %94 fit most images into 12KB/24KB block
                        end
                    end
                end
            end
        end
    end
    parsave(['image_sum_shift_',num2str(i_shift),'.mat'],image_sum,num) %save mean data of the image
end
toc