clear; %close all;clc
%demo script for PACBED_tilt(13.6mrad) CNN data prep
%Weizong Xu, July, 2017
%Please contact Prof. James LeBeau (jmlebeau@ncsu.edu) for simulation files.
%% get more images (generate signal images)
load('PACBED_STO_14mrad_data_1nm_resolution_mrad_tilt_38cat_corr3.mat','name_list','tilt_list','list_total')
img_num=1;%number of noisy image in the same condition
fpwd='D:\STEM_pattern_data\STO_PACBED_p227_AlexNet_t_cat38_4mrad';
num_center=5;%number of random images per condition
clear name_list
index_list=[0,10,25];%noise level from 0-25 (0-no noise, 25-noiest)
%%
for j=1:length(tilt_list)
    name_tilt{j}=['STO_',tilt_list{j},'_mrad'];
    if ~exist([fpwd,'\',name_tilt{j}],'dir')
        mkdir(fpwd,name_tilt{j})
        disp([fpwd,'\',name_tilt{j},' is created.'])
    end
end
%%
tic
for i_tilt=1:size(PACBED_data,1) %parfor compatible, but it is recommended to run multiple instances
    load('PACBED_STO_14mrad_data_1nm_resolution_mrad_tilt_38cat_corr3.mat','PACBED_data')
    PACBED_data=PACBED_data(i_tilt,:);%get all data from tilt
    num=uint32(0);
    image_sum=zeros(227,227,3);
    for i_list=5:length(list_total) %hard to get tilt angle thinner than 5nm
        for i_size=275:15:440 %set size range, value is the crop image pixel size
            PACBED_read=PACBED_data{1,i_list};
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
            image_size=size(PACBED_read,1);
            
            for i_inten=1.0 %set 1, i.e. no intensity scalling, training programm will do such job
                option.max_intensity_ratio=i_inten;
                for i_angle=0:2:88 % set rotation angle
                    disp([tilt_list{i_tilt},'_TiltCat_',num2str(i_list),'_Crop_',num2str(i_size),':',num2str(i_inten),' ', num2str(i_angle), ' #', num2str(num)])
                    PACBED_read_rot=imrotate(PACBED_read,i_angle,'bicubic','crop');
                    for i_index=1:length(index_list) %set noise level
                        option.noise_index=index_list(i_index);
                        for i_Center=1:num_center
                            img_center=[round(image_size/2), round(image_size/2)];
                            rand_x=image_size;rand_y=image_size;
                            while rand_x*rand_x+rand_y*rand_y>64 %set random shift range
                                rand_x=round((rand-0.5)*30);
                                rand_y=round((rand-0.5)*30);
                                option.signal_center(1)=rand_x+img_center(1);
                                option.signal_center(2)=rand_y+img_center(2);
                            end
                            option.distort_signal=1;
                            option.distort_y=(rand()-0.5)*0.02;%set random distorion via affine transformation
                            option.distort_y=(rand()-0.5)*0.02;
                            if i_Center==1 %be sure one training data is in the center
                                option.image_shift=[0,0];
                                option.signal_center=img_center;
                                option.distort_x=0;
                                option.distort_y=0;
                            end
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
                                filename=[fpwd,'/',name_tilt{i_tilt},'/Img_',num2str(num),'.jpg'];
                                if rand<0.5 %rand flip
                                    imwrite(uint8(flip(image_out{i_img},2)),filename, 'Quality',94) %94 fit most images into 12KB/24KB block
                                else
                                    imwrite(uint8(image_out{i_img}),filename, 'Quality',94)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    parsave(['image_sum_tilt_',num2str(i_tilt),'_cat38_14mrad.mat'],image_sum,num) %save mean data of the image
end
toc