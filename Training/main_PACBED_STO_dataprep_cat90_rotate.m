clear; %close all;clc
%demo script for PACBED_rotation CNN data prep
%Weizong Xu, July, 2017
%Please contact Prof. James LeBeau (jmlebeau@ncsu.edu) for simulation files.
%% get more images (generate signal images)
load('PACBED_STO_14mrad_data_1nm_resolution_corr3.mat')
PACBED_data1=PACBED_data;
load('PACBED_STO_19mrad_data_1nm_resolution_corr3.mat')
PACBED_data=[PACBED_data1; PACBED_data];
clear PACBED_data1
%%
img_num=1;
fpwd='C:\STEM_pattern_data\STO_PACBED_p227_rotate_13_19mrad';
num_center=10;
index_list=[0];
angle_list=-44:1:45;
for i=1:length(angle_list)
    if ~exist([fpwd,'\','STO_',num2str(angle_list(i)),'deg'],'dir')
        mkdir(fpwd,['STO_',num2str(angle_list(i)),'deg'])
        disp([fpwd,'\','STO_',num2str(angle_list(i)),'deg',' is created.'])
    end
end
%%
tic
parfor i_angle=-44:1:45
    image_sum=zeros(227,227,3);
    num=uint32(0);
    for i_list=1:length(list_total)
        for i_size=240:10:430 %for mbfit
            for i_cat=1:2
            PACBED_read=PACBED_data{i_cat,i_list};
            option=[];
            option.crop_size=i_size;%320
            rng('shuffle')
            option.rng_set=0; % run shuffle gloably to speed up
            option.noise_type='gaussian';%'poisson';
            option.rot_Angle=0;
            option.image_output_size=227;
            option.image_scale=255;
            option.noise_level=0;
            option.signal_level=10;
            option.chk_print=0;
            option.Gaus_blur_size=1;
            option.Avg_size=1;
            image_size=size(PACBED_read,1);
            
            for i_inten=1.0
                option.max_intensity_ratio=i_inten;
                disp([num2str(i_angle), 'deg-', name_list{i_list},':',num2str(i_inten),' ',num2str(i_size),' ', ' #', num2str(num)])
                PACBED_read_rot=imrotate(PACBED_read,i_angle,'bicubic','crop');
                for i_index=1:length(index_list)
                    option.noise_index=index_list(i_index);
                    for i_Center=1:num_center
                        img_center=[round(image_size/2), round(image_size/2)];
                        rand_x=image_size;rand_y=image_size;
                        while rand_x*rand_x+rand_y*rand_y>64
                            rand_x=round((rand-0.5)*30);
                            rand_y=round((rand-0.5)*30);
                            option.signal_center(1)=rand_x+img_center(1);
                            option.signal_center(2)=rand_y+img_center(2);
                        end
                        option.distort_signal=1;
                        option.distort_x=(rand()-0.5)*0.02;
                        option.distort_y=(rand()-0.5)*0.02;
                        if i_Center==1 %be sure one training data is in the center
                            option.image_shift=[0,0];
                            option.signal_center=img_center;
                            option.distort_x=0;
                            option.distort_y=0;
                        end
                        [ image_out, image_raw, ~ ] = PACBED_add_noise_gray( PACBED_read_rot, img_num, option );
                        
                        if i_index>0
                            img_num_eff=img_num;
                        else
                            img_num_eff=1;
                        end

                        for i_img=1:img_num_eff
                            num=num+1;
                            image_sum=image_sum+double(image_out{i_img});
                            filename=[fpwd,'/','STO_',num2str(i_angle),'deg','/Img_',num2str(num),'.jpg'];
                            imwrite(uint8(image_out{i_img}),filename, 'Quality',94)
                        end
                    end
                end
            end
            end
        end
    end
    parsave(['image_sum_rotate_',num2str(i_angle),'.mat'],image_sum,num)
end
toc