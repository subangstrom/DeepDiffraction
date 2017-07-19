function [ img_out, para_out, option ] = align_PACBED( data_list, net_str, option )
%align PACBED and crop into 227x227 size
%optional function, get rotational angle 
%Weizong, April, 

%initiate parameters
if ~exist('option','var')
    option.shift_enhance=0; %1-on, others off
    option.size_stable=0; %if all PACBEDs are in same size (2D map mode)
    option.rotation_stable=0; % output PACBEDs in same rotation angle (2D map mode)
    option.shift_avg=0; %average shift using 2D neighbour info
    option.crop_size_opt=330;
    option.output_rot_angle=nan; %reset output PACBED angle to value set, by default nan, i.e. not reset angle 
    option.mode='process_batch';
end

if isfield(option,'multiple_size_internal') % internal use, for fast thickness/tilt datafeed
    if option.multiple_size_internal==1
        [ img_out ] = reset_rotation_multiple_size_internal(data_list, option.center, option.crop_size, option.size_list, option.comp_size, option.rotation_angle, option.rotation_angle_out);
        return;
    end
end

%check version
curr_ver=version('-release');
if str2double(curr_ver(1:4))<2017
    warning('Matlab 2017a and later version is recommended.')
end

if ~isfield(option,'mode')
    option.mode='process_batch';
end

if strcmp(option.mode,'process_batch')
    [ img_out, para_out, option] = batch_process(data_list, net_str, option);
    return;
end

if strcmp(option.mode,'process_only_single')
    para_out.center=option.center;
    para_out.crop_size=option.crop_size;
    [ img_out, para_out.rot_angle] = reset_rotation(data_list, option.center, option.crop_size, option.rotation_angle, option.rotation_angle_out);
    return;
end


option.mode=lower(option.mode);
if strcmp(option.mode,'auto')
    option=[];
    option.rough_align=1;
    option.rotation=1;
    option.shift=1;
    option.size=1;
    option.crop_size_opt=320;
    option.accuracy_level=2; %0-fast, 1-normal, 2-accurate
    option.rotation_angle_out=nan;
    option.mode='auto';
end

if strcmp(option.mode(1:min(end,12)),'process_only')
    option.rough_align=0;
    option.rotation=0;
    option.shift=0;
    option.size=0;
    if ~isfield(option,'center') || ~isfield(option,'crop_size') || ~isfield(option,'rotation_angle')
        disp('Error! Not enough info provided (center/crop_size/rotation_angle). Return NULL.')
        img_out=[];
        para_out=[];
        return;
    end
end

if ~isfield(option, 'rotation_angle_out')
    option.rotation_angle_out=nan;
end

chk_rot_ini=0;
if ~isnan(option.rotation_angle_out)
    option.rotation=1;
    chk_rot_ini=1;
end

if ~isfield(option, 'rotation_angle')
    option.rotation_angle=[];
end

%load network
if ~isfield(option, 'rotation')
    option.rotation=0;
end
if option.rotation == 1
    try
        convnet_rotate=net_str.convnet_rotate;
    catch
        if chk_rot_ini==0
            disp('WARNING! No CNN net for rotation input. Disable rotation angle prediction.')
        end
        option.rotation=0;
    end
end
    
if ~isfield(option, 'shift')
    option.shift=0;
end
if option.shift == 1
    try
        convnet_shift=net_str.convnet_shift;
    catch
        disp('WARNING! No CNN net for shift input. Disable shift correction.')
        option.shift=0;
    end
end

if ~isfield(option, 'size')
    option.size=0;
end
if option.size == 1
    try
        convnet_size=net_str.convnet_size;
    catch
        disp('WARNING! No CNN net for size input. Disable size correction.')
        option.size=0;
    end
end


if ~isfield(option, 'crop_size_opt') || option.crop_size_opt<240 || option.crop_size_opt>400
    option.crop_size_opt=320;
else
    crop_size_opt=option.crop_size_opt;
end

if ~isfield(option, 'rough_align')
    option.rough_align=1; %enable rough align
end

if option.rough_align~=1 && (~isfield(option, 'center') || ~isfield(option, 'crop_size'))
    disp('Disable rough alignment, but not input center/crop_size. Enable rough alignment now.')
    option.rough_align=1;
end

if ~isfield(option, 'accuracy_level') || (option.accuracy_level~=0 && option.accuracy_level~=2)
    option.accuracy_level=1;
end

if option.accuracy_level==0
    iter_num_main=1;
    iter_shift=1;
    iter_size=1;
end

if option.accuracy_level==1
    iter_num_main=1;
    iter_shift=20;
    iter_size=10;
end

if option.accuracy_level==2
    iter_num_main=5;
    iter_shift=30;
    iter_size=20;
end

if ~iscell(data_list)
    data_list={data_list};
    chk_cell=0;
else
    chk_cell=1;
end

for i=1:length(data_list)
    if chk_cell==1
        disp(['Processing image #',num2str(i)])
    end
    %%load experiment data and rough crop
    PACBED_read=data_list{i};
    if ischar(PACBED_read)
        PACBED_read=imread(PACBED_read);
    end
    PACBED_read=single(PACBED_read(:,:,1));
    if max(max(PACBED_read))>256 %not uint8 style
        disp('Input data is set as 0-255 range, single precision.')
        PACBED_read=PACBED_read/max(max(PACBED_read))*255;
    end
%     PACBED_read=imresize(PACBED_read,0.5); %for testing
    %rough alignment
    if option.rough_align==1
        [ test_img, center, crop_size ] = rough_align_center( PACBED_read );
    else
        center=option.center;
        crop_size=option.crop_size;
        if iscell(center) 
            center_i=center{i};
        else
            center_i=center;
        end
        if iscell(crop_size)
            crop_size_i=crop_size{i};
        else
            crop_size_i=crop_size;
        end
        if ~strcmp(option.mode,'process_only')
            [ test_img ] = func_crop_image( PACBED_read, center_i, crop_size_i, 0 );
        end
    end
    
    for iter=1:iter_num_main
        %correct image shift
        if option.shift==1
            [ test_img, center ] = align_shift( test_img, PACBED_read, convnet_shift, center, crop_size, iter_shift );
        end
        %%correct image size
        if option.size==1
            [ test_img, crop_size ] = align_size( test_img, PACBED_read, convnet_size, center, crop_size, crop_size_opt, iter_size );
        end
    end
    
    %%analysis rotation angle
    if option.rotation == 1
        [ rot_angle ] = get_rotation_angle( test_img, convnet_rotate );
    else
        if iscell(option.rotation_angle)
            rot_angle=option.rotation_angle{i};
        else
            if isempty(option.rotation_angle)
                rot_angle=[];
            elseif length(option.rotation_angle)==1
                rot_angle=option.rotation_angle;
            else
                rot_angle=option.rotation_angle(i);
            end
        end   
    end

    %output figure if rotation angle for output is set
    if ~isnan(option.rotation_angle_out)
        if iscell(center) 
            center_i=center{i};
        else
            center_i=center;
        end
        if iscell(crop_size)
            crop_size_i=crop_size{i};
        else
            crop_size_i=crop_size;
        end
        [ test_img, rot_angle] = reset_rotation(PACBED_read, center_i, crop_size_i, rot_angle, option.rotation_angle_out);
    end


    %output data
    if chk_cell==1
        img_out{i}=test_img;
        if iscell(center) 
            para_out.center{i}=center{i};
        else
            para_out.center{i}=center;
        end
        if iscell(crop_size)
            para_out.crop_size{i}=crop_size{i};
        else
            para_out.crop_size{i}=crop_size;
        end
        para_out.rot_angle{i}=rot_angle;
    else
        img_out=test_img;
        para_out.center=center;
        para_out.crop_size=crop_size;
        para_out.rot_angle=rot_angle;
    end



    
end
end


%******************************FUNCTION*******************************
function [ img_out ] = func_crop_image( img_in, center, crop_size, show_fig )

try
    crop_img=img_in(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
catch
    %some padding added
    img_pad=zeros(size(img_in,1)+200, size(img_in,2)+200, size(img_in,3),'uint8');
    img_pad(101:end-100,101:end-100,:)=img_in;
    center=center+100;
    crop_img=img_pad(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
end
    img_out=imresize(crop_img,[227 227]);

if show_fig==1
    figure;imshow(img_out(:,:,1),[])
end

end

function [ test_img, center, crop_size ] = rough_align_center( PACBED_read )
PACBED_read(PACBED_read==0)=1; %avoid nan
sum_x=sum(PACBED_read,2);
sum_y=sum(PACBED_read,1);
% figure;plot(sum_x)
% figure;plot(sum_y)
% [~,pos_x]=max(sum_x);
% [~,pos_y]=max(sum_y);
para_gauss_x = gaussianFit([1:size(PACBED_read,2)]',sum_x);
para_gauss_y = gaussianFit([1:size(PACBED_read,1)]',sum_y);

crop_size=(para_gauss_x(2)+para_gauss_y(2))*1.0;
% center=[pos_x, pos_y];
center=round([para_gauss_x(1),para_gauss_y(1)]);
% center0=center;
[ test_img ] = func_crop_image( PACBED_read, center, crop_size, 0 );
test_img=uint8(cat(3,test_img,test_img,test_img));
% [ test_img11 ] = func_crop_image( test_img, [114,114], 150, 1 );
end

function [ test_img, center ] = align_shift( test_img, PACBED_read, convnet_shift, center, crop_size, iter_shift )
%%search for center
iter_num=0;
while iter_num<iter_shift
    shift_x = floor((classify_PACBED_shift( convnet_shift, test_img, 0)-classify_PACBED_shift( convnet_shift, test_img, 180))/2);
    shift_y = floor((classify_PACBED_shift( convnet_shift, test_img, 90)-classify_PACBED_shift( convnet_shift, test_img, -90))/2);
    if abs(shift_x)+abs(shift_y)>0
        center=center+[-shift_x, +shift_y];%tested for correctness, because coordinate defined in image is different
        [ test_img ] = func_crop_image( PACBED_read, center, crop_size, 0 );
        test_img=uint8(cat(3,test_img,test_img,test_img));
    else
        break;
    end
    iter_num=iter_num+1;
end
%%quick check
% figure;imshow(test_img(:,:,1),[])
% func_crop_image( test_img, [114,114], 150, 1 );

end

function [ test_img, crop_size ] = align_size( test_img, PACBED_read, convnet_size, center, crop_size, crop_size_opt, iter_size )
iter_num=0;
while iter_num<iter_size
    size_PACBED = classify_PACBED_size( convnet_size, test_img);
    if size_PACBED~=crop_size_opt
        crop_size=crop_size*crop_size_opt/size_PACBED;
        [ test_img ] = func_crop_image( PACBED_read, center, crop_size, 0 );
        test_img=uint8(cat(3,test_img,test_img,test_img));
    else
        break;
    end
    iter_num=iter_num+1;
end
end

function [ shift_x ] = classify_PACBED_shift( convnet_shift, test_img, rot_angle )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('rot_angle','var') || rot_angle==0
    Ypred_x = classify(convnet_shift,test_img);
else
    Ypred_x=classify(convnet_shift,imrotate(test_img,round(rot_angle/90)*90));
end

shift_x=char(Ypred_x);
shift_x=shift_x(11:end);
shift_x=str2double(shift_x);

end

function [ size ] = classify_PACBED_size( convnet_size, test_img, rot_angle )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('rot_angle','var') || rot_angle==0
    Ypred_x = classify(convnet_size,test_img);
else
    Ypred_x=classify(convnet_size,imrotate(test_img,round(rot_angle/90)*90));
end

size=char(Ypred_x);
size=size(10:end);
size=str2double(size);

end

function [ rot_angle ] = get_rotation_angle( test_img, convnet_rotate )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Ypred_x = classify(convnet_rotate,test_img);

rot_angle=char(Ypred_x);
rot_angle=rot_angle(5:end-3);
rot_angle=-str2double(rot_angle); %define positive - clockwise

end

function [ test_img, rotation_angle_out] = reset_rotation(img, center, crop_size, rot_angle, rotation_angle_out)
%Crop from rotated image with given angle
%define positive angle - clockwise
%Weizong Xu, April 2017.

img_center=round((size(img)-1)/2);
center_diff=center-img_center;
theta=rot_angle-rotation_angle_out;
tmp_cos=cosd(theta);
tmp_sin=sind(theta);
rot_M=[tmp_cos -tmp_sin; tmp_sin tmp_cos];
center_diff_rot=rot_M*center_diff';  
img_rot=imrotate(img,theta,'crop');
img_center_rot=round(size(img_rot)-1)/2;
center_rot=round(center_diff_rot'+img_center_rot);
[ test_img ] = func_crop_image( img_rot, center_rot, crop_size, 0 );
if isa(test_img,'uint8')
    test_img=uint8(cat(3,test_img,test_img,test_img));
end

end

function [ test_img ] = reset_rotation_multiple_size_internal(img, center, crop_size, size_list, comp_size, rot_angle, rotation_angle_out)
%only for internal use
img_center=round((size(img)-1)/2);
center_diff=center-img_center;
theta=rot_angle-rotation_angle_out;
tmp_cos=cosd(theta);
tmp_sin=sind(theta);
rot_M=[tmp_cos -tmp_sin; tmp_sin tmp_cos];
center_diff_rot=rot_M*center_diff';  
img_rot=imrotate(img,theta,'crop');
img_center_rot=round(size(img_rot)-1)/2;
center_rot=round(center_diff_rot'+img_center_rot);
test_img=cell(length(size_list),length(comp_size));
for i=1:length(size_list)
    for j=1:length(comp_size)
        crop_size_mod=crop_size*(1+(size_list(i)+comp_size(j))/100);
        test_img{i,j}=func_crop_image( img_rot, center_rot, crop_size_mod, 0 );
    end
end

if isa(test_img{1,1},'uint8')
    for i=1:length(size_list)
        for j=1:length(comp_size)
            test_img{i,j}=uint8(cat(3,test_img{i,j},test_img{i,j},test_img{i,j}));
        end
    end
end

end


function [ database_crop, para_out, option] = batch_process(database_in, net_str, option)
% tic
%Setup parameters
%setup database
size_input=size(database_in);
if length(size_input)>2
    disp('Database in 3 or more dimension is not supported. Please reshape to 1D dimension to input.')
    option.error=1;
    database_crop=[];
    para_out=[];
    return;
end

if ~isfield(net_str,'convnet_size') || ~isfield(net_str,'convnet_shift')
    msgbox('Alignment failed! No size or shift network input!')
    error('Alignment failed! No size or shift network input!')
end

if ~isfield(net_str,'convnet_rotate')
    warning('No rotation network input! Set rotation angle output as 0.')
    net_str.convnet_rotate=[];
end

if size_input(2)>1 && size_input(1)>1
    database=reshape(database_in',[size_input(1)*size_input(2) 1]);
    chk_reshape=1;
else
    database=database_in;
    chk_reshape=0;
end


if ~isfield(option, 'shift_enhance')
    shift_enhance=1; %1-on, others off
else
    shift_enhance=option.shift_enhance;
end

if ~isfield(option, 'crop_size_opt')
    crop_size_opt=330;
else
    crop_size_opt=option.crop_size_opt;
end

if ~isfield(option, 'size_stable')
    if chk_reshape==1
        size_stable=1; %if all PACBEDs are in same size (2D map mode)
    else
        size_stable=0; %1D input not stablized the pattern size
    end
else
    size_stable=option.size_stable;
end

if ~isfield(option, 'rotation_stable')
    if chk_reshape==1
        rotation_stable=1; % output PACBEDs in same rotation angle (2D map mode)
    else
        rotation_stable=0; % 1D input not stablized the pattern rotation angle
    end
else
    rotation_stable=option.rotation_stable;
end

if ~isfield(option, 'iter_num_max')
    iter_num_max=6; %max iter for shift
else
    iter_num_max=round(option.iter_num_max);
end

if ~isfield(option, 'iter_size_stable_max')
    iter_size_stable_max=3; %max iter for size
else
    iter_size_stable_max=round(option.iter_size_stable_max);
end

if ~isfield(option, 'Int_cutoff')
    Int_cutoff=uint8(50); %reduce excess BG contrast
else
    Int_cutoff=uint8(option.Int_cutoff);
end

if ~isfield(option, 'database_crop_ratio')
    database_crop_ratio=4; %reduce size for gaussian fit
else
    database_crop_ratio=option.database_crop_ratio;
end

if ~isfield(option, 'shift_avg')
    shift_avg=0; %avg size for shift
else
    shift_avg=round(option.shift_avg);
end

if ~isfield(option, 'output_rot_angle')
    output_rot_angle=nan; %nan - no reset angle
else
    output_rot_angle=option.output_rot_angle;
end

if ~isfield(option, 'cal_gpu')
    cal_gpu=1; %enable gpu
else
    cal_gpu=option.cal_gpu;
end

if cal_gpu==1
    cal_gpu=GPU_check(net_str.convnet_size)+GPU_check(net_str.convnet_shift)-1;
    if cal_gpu~=1
        msgbox('GPU is disabled, use CPU to run CNN. It could take much longer time to finish.');
    end
end

if cal_gpu~=1
    disp('GPU is disabled, use CPU to run CNN. It could take much longer time to finish.');
end

%rough_align_center(database);
database_crop=database;
database_pred=zeros([net_str.convnet_shift.Layers(1,1).InputSize,length(database)*2],'uint8');
img_size=net_str.convnet_shift.Layers(1,1).InputSize(1);
database_center=cell(length(database),1);
database_size=zeros(length(database),1);
Int_cutoff=uint8(Int_cutoff);
% tic
for i=1:length(database)
    if ~isa(database{i},'uint8')
        database{i}=double(database{i});
        database{i}=uint8(database{i}/max(max(max(database{i})))*255);
    end
    if size(database{i},1)<=300 || size(database{i},2)<=300
        center=round([size(database{i},1)/2, size(database{i},2)/2]);
        crop_size=min([size(database{i},1), size(database{i},2)])-2;
    else
        PACBED_read=imresize(database{i}(:,:,1),1/database_crop_ratio);
        PACBED_read=double(PACBED_read);
        PACBED_read(PACBED_read==0)=1; %avoid nan
        sum_x=sum(PACBED_read,2);
        sum_y=sum(PACBED_read,1);
    %     [~,pos_x]=max(sum_x);
    %     [~,pos_y]=max(sum_y);
        para_gauss_x = gaussianFit((1:size(PACBED_read,2))',sum_x);
        para_gauss_y = gaussianFit((1:size(PACBED_read,1))',sum_y);
        crop_size=(para_gauss_x(2)+para_gauss_y(2))*database_crop_ratio*1.25;
        center=round([para_gauss_x(1),para_gauss_y(1)])*database_crop_ratio;
    %     center=[pos_x, pos_y]*database_crop_ratio;
    end
    database_center{i}=center;
    database_size(i)=crop_size;
    crop_img=database{i}(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
    database_crop{i}=imresize(crop_img,[img_size img_size]);
    tmp_img=cat(3,database_crop{i},database_crop{i},database_crop{i})-Int_cutoff;
    database_pred(:,:,:,i)=tmp_img;
    database_pred(:,:,:,i+length(database))=rot90(tmp_img);
    if shift_enhance==1
        database_pred(:,:,:,i+length(database)*2)=rot90(tmp_img,2);
        database_pred(:,:,:,i+length(database)*3)=rot90(tmp_img,3);
    end
end
% toc
% Img_series_Looper(database_crop,name_list_1D);

Class_table_shift=net_str.convnet_shift.Layers(end,1).ClassNames;
Score_table_shift=zeros(length(Class_table_shift),1);
for i=1:length(Class_table_shift)
    Score_table_shift(i,1)=str2double(Class_table_shift{i,1}(11:end));
end
Class_table_size=net_str.convnet_size.Layers(end,1).ClassNames;
Score_table_size=zeros(length(Class_table_size),1);
for i=1:length(Class_table_size)
    Score_table_size(i,1)=str2double(Class_table_size{i,1}(10:end));
end

if ~isempty(net_str.convnet_rotate)
    Class_table_rotate=net_str.convnet_rotate.Layers(end,1).ClassNames;
    Score_table_rotate=zeros(length(Class_table_rotate),1);
    for i=1:length(Class_table_rotate)
        Score_table_rotate(i,1)=str2double(Class_table_rotate{i,1}(5:end-3));
    end
else
    Score_table_rotate=nan;
end

for iter_num=1:iter_num_max
    if cal_gpu==1
        [~,scores_shift] = classify(net_str.convnet_shift,database_pred);
    else
        [~,scores_shift] = classify(net_str.convnet_shift,database_pred,'ExecutionEnvironment','cpu');
    end
    [~,scores_shx]=max(scores_shift(1:length(database),:),[],2); %use most probable value
    scores_shx(:)=Score_table_shift(scores_shx(:));
    [~,scores_shy]=max(scores_shift(length(database)+1:length(database)*2,:),[],2);
    scores_shy(:)=Score_table_shift(scores_shy(:));
    if shift_enhance==1
        [~,scores_shx_inv]=max(scores_shift(length(database)*2+1:length(database)*3,:),[],2); %use most probable value
        scores_shx_inv(:)=Score_table_shift(scores_shx_inv(:));
        scores_shx=round((scores_shx-scores_shx_inv)*0.5);
        [~,scores_shy_inv]=max(scores_shift(length(database)*3+1:length(database)*4,:),[],2); %use most probable value
        scores_shy_inv(:)=Score_table_shift(scores_shy_inv(:));
        scores_shy=round((scores_shy-scores_shy_inv)*0.5);
    end
    if iter_num==1
        for i=1:length(database)
            database_center{i}=database_center{i}+[-scores_shx(i), +scores_shy(i)];
            center=database_center{i};
            crop_size=database_size(i);
            try
                crop_img=database{i}(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
            catch
                %some padding added
                img_pad=zeros(size(database{i},1)+200, size(database{i},2)+200, size(database{i},3),'uint8');
                img_pad(101:end-100,101:end-100,:)=database{i};
                center=center+100;
                crop_img=img_pad(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
            end
            
            database_crop{i}=imresize(crop_img,[img_size img_size]);
            tmp_img=cat(3,database_crop{i},database_crop{i},database_crop{i})-Int_cutoff;
            database_pred(:,:,:,i)=tmp_img;
            database_pred(:,:,:,i+length(database))=rot90(tmp_img);
            if shift_enhance==1
                database_pred(:,:,:,i+length(database)*2)=rot90(tmp_img,2);
                database_pred(:,:,:,i+length(database)*3)=rot90(tmp_img,3);
            end
        end
        continue;
    end
    if size_stable==1 && iter_num<=iter_size_stable_max %run at iter_num==2,3
        if cal_gpu==1
            [~,scores_size] = classify(net_str.convnet_size,database_pred(:,:,:,1:length(database)));
        else
            [~,scores_size] = classify(net_str.convnet_size,database_pred(:,:,:,1:length(database)),'ExecutionEnvironment','cpu');
        end
        [~,scores_si]=max(scores_size,[],2);
        scores_si(:)=Score_table_size(scores_si(:));
    end
    if size_stable==1 && iter_num>iter_size_stable_max
        database_size(:)=median(database_size);
        scores_si(:)=crop_size_opt;%median(scores_si);
    end
    if size_stable~=1
        if cal_gpu==1
            [~,scores_size] = classify(net_str.convnet_size,database_pred(:,:,:,1:length(database)));
        else
            [~,scores_size] = classify(net_str.convnet_size,database_pred(:,:,:,1:length(database)),'ExecutionEnvironment','cpu');
        end
        [~,scores_si]=max(scores_size,[],2);
        scores_si(:)=Score_table_size(scores_si(:));
    end
    for i=1:length(database)
        database_center{i}=database_center{i}+[-scores_shx(i), +scores_shy(i)];
        center=database_center{i};
        database_size(i)=round(database_size(i)*crop_size_opt/scores_si(i));
        crop_size=database_size(i);
        try
            crop_img=database{i}(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
        catch
            %some padding added
            img_pad=zeros(size(database{i},1)+200, size(database{i},2)+200, size(database{i},3),'uint8');
            img_pad(101:end-100,101:end-100,:)=database{i};
            center=center+100;
            crop_img=img_pad(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
        end
        database_crop{i}=imresize(crop_img,[img_size img_size]);
        tmp_img=cat(3,database_crop{i},database_crop{i},database_crop{i})-Int_cutoff;
        database_pred(:,:,:,i)=tmp_img;
        database_pred(:,:,:,i+length(database))=rot90(tmp_img);
        if shift_enhance==1
            database_pred(:,:,:,i+length(database)*2)=rot90(tmp_img,2);
            database_pred(:,:,:,i+length(database)*3)=rot90(tmp_img,3);
        end
    end

    chk_si=max(abs(scores_si-median(scores_si)));
    chk_sh=max([scores_shx;scores_shy]);
    if chk_si==0 && chk_sh<=1
        break;
    end
end

if ~isempty(net_str.convnet_rotate)
    if cal_gpu==1
        [~,scores_rotate] = classify(net_str.convnet_rotate,database_pred(:,:,:,1:length(database)));
    else
        [~,scores_rotate] = classify(net_str.convnet_rotate,database_pred(:,:,:,1:length(database)),'ExecutionEnvironment','cpu');
    end
    [~,scores_r]=max(scores_rotate,[],2);
    scores_r(:)=-Score_table_rotate(scores_r(:)); %positive at clockwise rot
else
    scores_r=zeros(length(database),1);
end

% toc
%Generate aligned database and check
if chk_reshape==1 && shift_avg>1
    center_2D=reshape(database_center,[size_input(2) size_input(1)])';
    shiftx_2D=zeros(size_input(1),size_input(2));
    shifty_2D=shiftx_2D;
    for ix=1:size_input(1)
        for jy=1:size_input(2)
            shiftx_2D(ix,jy)=center_2D{ix,jy}(1);
            shifty_2D(ix,jy)=center_2D{ix,jy}(2);
        end
    end
    shiftx_2D=round(avg_image(shiftx_2D,shift_avg));
    shifty_2D=round(avg_image(shifty_2D,shift_avg));
    shiftx_1D=reshape(shiftx_2D',[size_input(1)*size_input(2) 1]);
    shifty_1D=reshape(shifty_2D',[size_input(1)*size_input(2) 1]);
    for ix=1:size_input(1)*size_input(2)
        database_center{ix}=[shiftx_1D(ix),shifty_1D(ix)];
    end
    for i=1:length(database)
        center=database_center{i};
        crop_size=database_size(i);
        crop_img=database{i}(center(1)-floor(crop_size/2)+1:center(1)+floor(crop_size/2),center(2)-floor(crop_size/2)+1:center(2)+floor(crop_size/2),:);
        database_crop{i}=imresize(crop_img,[img_size img_size]);
    end
end

scores_r_mid=median(scores_r);
if rotation_stable==1
    scores_r(:)=scores_r_mid;
end

if ~isnan(output_rot_angle)
    for i=1:length(database)
        [ database_crop{i}, ~] = reset_rotation(database{i}, database_center{i}, database_size(i), scores_r(i), output_rot_angle);
    end
end

%Output aligned parameter
para_out=cell(length(database),1);
for i=1:length(database)
    para_out{i,1}.center=database_center{i};
    para_out{i,1}.crop_size=database_size(i);
    para_out{i,1}.rot_angle=scores_r(i);
end
% toc
if chk_reshape==1
    database_crop=reshape(database_crop,[size_input(2) size_input(1)])';
    para_out=reshape(para_out,[size_input(2) size_input(1) ])';
end

end