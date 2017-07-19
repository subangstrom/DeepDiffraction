function [ image_out, image_raw, option ] = PACBED_add_noise_gray( PACBED_read, img_num, option )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% chk_gray2rgb=0;
if size(PACBED_read,3)>1
    PACBED_read=rgb2gray(PACBED_read);
%     chk_gray2rgb=1;
end

if ~exist('img_num','var') || img_num<1
    img_num=1;
end
img_num=floor(img_num);

if ~exist('option','var')
    option.autoGen=1; %auto generation;
end
    
if ~isfield(option,'rng_set')
    option.rng_set=1; %use shuffle in rng
end

if ~isfield(option,'image_scale')
    option.image_scale=255; %for 255
end

if ~isfield(option,'BG_level')
    option.BG_level=0; %background level 
end

if ~isfield(option,'noise_level')
    option.noise_level=0; %noise level 
end

if ~isfield(option,'noise_index')
    option.noise_index=0; %noise index from 0-25
end

if ~isfield(option,'signal_level')
    option.signal_level=10; %signal level
end
%calculation take signal+noise+background and noise index into consideration

if ~isfield(option,'max_intensity_ratio') || option.max_intensity_ratio<0 || option.max_intensity_ratio>1
    option.max_intensity_ratio=1.0; % 100%
end

if ~isfield(option,'signal_type')
    option.signal_type='Gaussian'; %produce 2D Gaussian signal
end

if ~isfield(option,'noise_type') % Gaussian or Poisson
    option.noise_type='Gaussian'; %produce Gaussian distributed noise
end

if ~isfield(option,'rot_Angle')
    option.rot_Angle=0;
end
if option.rot_Angle~=0
    PACBED_read=imrotate(PACBED_read,option.rot_Angle,'bicubic','crop');
end

if ~isfield(option,'signal_center')
    option.signal_center=[round(size(PACBED_read,1)/2), round(size(PACBED_read,2)/2)]; %Gaussian Sigma=25;
end

if ~isfield(option,'image_shift')
    option.image_shift=[0, 0];
end

if ~isfield(option,'Gaus_blur_size')
    option.Gaus_blur_size=1; %1-no Gaussian Filter
end
option.Gaus_blur_size=floor(option.Gaus_blur_size);

if ~isfield(option,'Avg_size')
    option.Avg_size=1; %1 - no Average Filter
end
option.Avg_size=floor(option.Avg_size);

if ~isfield(option,'chk_print')
    option.chk_print=0; %do not display figure
end

if ~isfield(option,'distort_signal')
    option.distort_signal=0; %do not distort image
end

if ~isfield(option,'distort_x')
    option.distort_x=0;
end

if ~isfield(option,'distort_y')
    option.distort_y=0;
end

if option.distort_signal~=0 && option.distort_x==0 && option.distort_y==0
    option.distort_x=(rand()-0.5)*0.002;
    option.distort_y=(rand()-0.5)*0.002;
end

if option.distort_signal~=0
    Am=[1 option.distort_y 0; option.distort_x 1 0; 0 0 1];
    tform = affine2d(Am);
    A_z=PACBED_read*0;
    A_signal_tmp = imwarp(PACBED_read,tform);
    ims=floor(size(PACBED_read)/2);
    ims_diff=size(PACBED_read)-ims*2;
    tmp=floor(size(A_signal_tmp)/2);
    A_z(ims(1)-min(ims(1),tmp(1))+1:ims(1)+min(ims(1),tmp(1))+ims_diff(1),ims(2)-min(ims(2),tmp(2))+1:ims(2)+min(ims(2),tmp(2))+ims_diff(2))=...
    A_signal_tmp(tmp(1)-min(ims(2),tmp(1))+1:tmp(1)+min(ims(2),tmp(1))+ims_diff(1),tmp(2)-min(ims(2),tmp(2))+1:tmp(2)+min(ims(2),tmp(2))+ims_diff(2));
    PACBED_read=A_z;
end

if ~isfield(option,'crop_size') || option.crop_size<=0
    option.crop_size=min([320,size(PACBED_read,1),size(PACBED_read,2)]);
end
if option.crop_size<size(PACBED_read,1)
    crop_center=option.signal_center;%round(size(PACBED_read)/2);
    PACBED_read=PACBED_read(crop_center(1)-round(option.crop_size/2)+1:crop_center(1)+round(option.crop_size/2),crop_center(2)-round(option.crop_size/2)+1:crop_center(2)+round(option.crop_size/2),:);
end

if ~isfield(option,'image_output_size') || option.image_output_size<=0
    option.image_output_size=227;
end
image_size=option.image_output_size;
if image_size~=size(PACBED_read,1) || image_size~=size(PACBED_read,2)
    image_raw=imresize(PACBED_read,[image_size image_size]);
    PACBED_resize=double(image_raw);
else
    image_raw=PACBED_read;
    PACBED_resize=double(PACBED_read);
end

if option.image_shift(1)~=0 && option.image_shift(2)~=0
    image_shift=option.image_shift;
    img_tmp=PACBED_resize;
    img_shift=img_tmp*0;
    img_shift(max(image_shift(1)+1,1):min(end,end+image_shift(1)),max(image_shift(2)+1,1):min(end,end+image_shift(2)))=...
        img_tmp(max(1,-image_shift(1)+1):min(end-image_shift(1),end),max(1,-image_shift(2)+1):min(end-image_shift(2),end));
    PACBED_resize=img_shift;
end

if option.rng_set==1
    rng('shuffle')
end

amp_noise=option.image_scale*option.max_intensity_ratio*(option.noise_level/(option.noise_level+option.signal_level+option.BG_level));
amp_BG=option.image_scale*option.max_intensity_ratio*(option.BG_level/(option.noise_level+option.signal_level+option.BG_level));
amp_signal=option.image_scale*option.max_intensity_ratio*(option.signal_level/(option.noise_level+option.signal_level+option.BG_level));
image_out=cell(img_num,1);

for i=1:img_num
    A_rand=rand(image_size,image_size)*amp_noise;
    A_base=ones(image_size,image_size)*amp_BG;
    A_signal=PACBED_resize(:,:,1)*(amp_signal/max(max(max(PACBED_resize))));
    A_tot=uint8(A_rand+A_base+A_signal);
    
    tot_integ_num=max(25-option.noise_index+1,1);
    if tot_integ_num>=1 && option.noise_index>0
        I_noise_tot=double(imnoise(uint8(A_tot),option.noise_type))/tot_integ_num;
        for integ=1:tot_integ_num-1
            I_noise_tot=I_noise_tot+double(imnoise(uint8(A_tot),option.noise_type))/tot_integ_num;
        end
    else
        I_noise_tot=A_tot;
    end
    image_out{i}=uint8(I_noise_tot);
    
    if option.chk_print==1 && i<=5
        figure;imagesc(image_out{i},[0 option.image_scale]);axis image;colormap(gray)
    end
    
    if (option.Gaus_blur_size>=2)
        image_out{i}=imgaussfilt(image_out{i},option.Gaus_blur_size);
        if option.chk_print==1 && i<=5
            figure;imagesc(image_out{i},[0 option.image_scale]);axis image;colormap(gray)
        end
    end
    
    if (option.Avg_size>=2)
        image_out{i}=avg_image(image_out{i},option.Avg_size);
        if option.chk_print==1 && i<=5
            figure;imagesc(image_out{i},[0 option.image_scale]);axis image;colormap(gray)
        end
    end
end

% if chk_gray2rgb==1 %need convert back to rgb
%     for i=1:img_num
%         image_out{i}=cat(3,image_out{i},image_out{i},image_out{i});
%     end
%     image_raw=cat(3,image_raw,image_raw,image_raw);
% end

end

