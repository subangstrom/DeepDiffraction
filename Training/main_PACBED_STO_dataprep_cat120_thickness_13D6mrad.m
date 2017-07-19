clear; %close all;clc
%demo script for PACBED_thickness(13.6mrad) CNN data prep
%Weizong Xu, July, 2017
%% get more images (generate signal images)
load('PACBED_STO_14mrad_data_1nm_resolution_mrad_tilt_38cat_corr3.mat')
img_num=1;
fpwd='E:\STEM_pattern_data\STO_PACBED_p227_14mrad_thickness';
num_center=5;
name_list_all=name_list;
clear name_list
index_list=[0,10,25];
%%
for i=1:length(list_total)
    name_list{i}=['STO_',num2str(i),'nm'];
    if ~exist([fpwd,'\',name_list{i}],'dir')
        mkdir(fpwd,name_list{i})
        disp([fpwd,'\',name_list{i},' is created.'])
    end
end
%%
tic
for i_list=1:120%120
    num=uint32(0);
    image_sum=zeros(227,227,3);
    for i_tilt=1:size(PACBED_data,1)
        for i_size=245:15:410
            PACBED_read=PACBED_data{i_tilt,i_list};
            option=[];
            option.crop_size=i_size;
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
                for i_angle=0:2:88
                    disp([name_list{i_list},'_TiltCat_',num2str(i_tilt),'_Crop_',num2str(i_size),':',num2str(i_inten),' ', num2str(i_angle), ' #', num2str(num)])
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
                            [ image_out, ~, ~ ] = PACBED_add_noise_gray( PACBED_read_rot, img_num, option );
                            
                            if i_index>0
                                img_num_eff=img_num;
                            else
                                img_num_eff=1;
                            end
                            
                            for i_img=1:img_num_eff
                                num=num+1;
                                image_sum=image_sum+double(image_out{i_img});
                                filename=[fpwd,'/',name_list{i_list},'/Img_',num2str(num),'.jpg'];
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
    parsave(['image_sum_thickness_',num2str(i_list),'_14mrad.mat'],image_sum,num)
end
toc
exit