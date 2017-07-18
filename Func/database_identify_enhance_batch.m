function [ thickness_determine, tilt_determine, data_out_all, tilt_align_cell, option] = database_identify_enhance_batch( database, net_str, para_in, option )
%Identify PACBED thickness and tilt using enhanced search
%Very fast speed using parallel batch
%Weizong Xu, May, 2017

%check version
curr_ver=version('-release');
if str2double(curr_ver(1:4))<2017
    warning('Matlab 2017a and later version is recommended.')
end

%initiate parameters
if ~exist('option','var')
    option.thickness_comp_ratio=0.70;
    option.angle_list=0:30:90;%0:30:90;
    option.size_list=-5:5:5;%-5:5:0;
    option.cal_gpu=1; % 1-gpu, 0-cpu
    option.thickness_stable=0; %stable thickness with +- input range from 2D data 
    option.comp_size_thickness=-5; %size compensation (%) in thickness calculation
    option.comp_size_tilt=5; %size compensation (%) in tilt calculation
    option.symmetry=4; %STO [001] four-fold; STO [011] two-fold
    option.aperture=[];
end

if ~isfield(option, 'thickness_comp_ratio')
    thickness_comp_ratio=0.70;
else
    thickness_comp_ratio=option.thickness_comp_ratio;
end

if ~isfield(option, 'angle_list')
    angle_list=0:30:90;
else
    angle_list=option.angle_list;
end

if ~isfield(option, 'size_list')
    size_list=-5:5:5;
else
    size_list=option.size_list;
end

if ~isfield(option, 'cal_gpu')
    cal_gpu=1; %enable gpu
else
    cal_gpu=option.cal_gpu;
end

if ~isfield(option, 'comp_size_thickness')
    comp_size_thickness=-5; %size compensation for optimal size in thickness calculation
else
    comp_size_thickness=option.comp_size_thickness;
end

if ~isfield(option, 'comp_size_tilt')
    comp_size_tilt=5; %size compensation for optimal size in tilt calculation
else
    comp_size_tilt=option.comp_size_tilt;
end

if cal_gpu==1
    cal_gpu=GPU_check(net_str.convnet_thickness)+GPU_check(net_str.convnet_tilt)-1;
end

if cal_gpu~=1
    disp('GPU is not working, use CPU to run CNN. It could take much longer time to finish.');
    msgbox('GPU is not working, use CPU to run CNN. It could take much longer time to finish.');
end

if ~isfield(option, 'thickness_stable')
    thickness_stable=0;
else
    thickness_stable=option.thickness_stable;
end

if ~isfield(option, 'aperture')
    cal_mask=[];
else
    cal_mask=option.aperture;
    mask_width=max(cal_mask)-min(cal_mask);
end

if ~isfield(option, 'symmetry')
    tilt_symmetry=4; %STO [001] four-fold; STO [011] two-fold
else
    tilt_symmetry=option.symmetry;
    if tilt_symmetry~=4 && tilt_symmetry~=2
        disp('The symmetry setting is currently not supported. Auto change it to 2-fold symmetry to continue the processing.')
        tilt_symmetry=2;
    end
end

if ~isfield(option,'symmetry_rot_offset')
    symmetry_rot_offset=0;
else
    symmetry_rot_offset=option.symmetry_rot_offset;
end

size_input=size(database);
if length(size_input)>2
    disp('Database in 3 or more dimension is not supported. Please reshape to 1D dimension to input.')
    thickness_determine=[];
    tilt_determine=[];
    data_out_all=[];
    tilt_align_cell=[];
    return;
end

size_check=size_input-size(para_in);
if size_check(1)+size_check(2)>0
    disp('The dimension of input image database and calibration are not consistent. Return NULL.')
    thickness_determine=[];
    tilt_determine=[];
    data_out_all=[];
    tilt_align_cell=[];
    return;
end

if size_input(2)>1 && size_input(1)>1
    database=reshape(database',[size_input(1)*size_input(2) 1]);
    para_in=reshape(para_in',[size_input(1)*size_input(2) 1]);
    chk_reshape=1;
else
    chk_reshape=0;
end

%data prepare for NN prediction
num=0;
total_img_count=length(database);
for i=1:length(database)
    if ~isa(database{i},'uint8')
        database{i}=double(database{i});
        database{i}=uint8(database{i}/max(max(max(database{i})))*255);
    end
end
big_database=zeros(227,227,3,total_img_count*length(angle_list)*length(size_list),'uint8');
big_database_tilt=big_database;
database_tilt_align=zeros(227,227,3,total_img_count,'uint8');
option1.multiple_size_internal=1;
if isempty(cal_mask)
    if comp_size_thickness==comp_size_tilt
        for ni=1:total_img_count
            option1.center=para_in{ni}.center;
            option1.rotation_angle=para_in{ni}.rot_angle;
            for i=angle_list
                option1.size_list=size_list;
                option1.comp_size=comp_size_thickness;
                option1.crop_size=para_in{ni}.crop_size;
                option1.rotation_angle_out=i;
                [ img_align1 ] = align_PACBED( database{ni}, [], option1 );
                for j=1:length(size_list)
                    num=num+1;
                    big_database(:,:,:,num)=img_align1{j,1};
                    if i==0 && size_list(j)==0
                        database_tilt_align(:,:,:,ni)=img_align1{j,1};
                    end
                end
            end
        end
        big_database_tilt=big_database;
    else
        for ni=1:total_img_count
            option1.center=para_in{ni}.center;
            option1.rotation_angle=para_in{ni}.rot_angle;
            for i=angle_list
                option1.size_list=size_list;
                option1.comp_size=[comp_size_thickness,comp_size_tilt];
                option1.crop_size=para_in{ni}.crop_size;
                option1.rotation_angle_out=i;
                [ img_align1 ] = align_PACBED( database{ni}, [], option1 );
                for j=1:length(size_list)
                    num=num+1;
                    big_database(:,:,:,num)=img_align1{j,1};
                    big_database_tilt(:,:,:,num)=img_align1{j,2};
                    if i==0 && size_list(j)==0
                        database_tilt_align(:,:,:,ni)=img_align1{j,2};
                    end
                end
            end
        end
    end

    if isempty(angle_list(angle_list==0)) || isempty(size_list(size_list==0)) %input not contain 0 degree/tilt
        for ni=1:total_img_count
            option11.mode='process_only_single';
            option11.center=para_in{ni}.center;
            option11.rotation_angle=para_in{ni}.rot_angle;
            option11.crop_size=para_in{ni}.crop_size*(1+comp_size_thickness/100);
            option11.rotation_angle_out=0;
            [ img_align2 ] = align_PACBED( database{ni}, [], option11 );
            database_tilt_align(:,:,:,ni)=img_align2;
        end
    end
else
    [rm_matrix, rm_matrix_inv] = mask_ring( ones(227,227), cal_mask, nan );
    if comp_size_thickness==comp_size_tilt
        for ni=1:total_img_count
            option1.center=para_in{ni}.center;
            option1.rotation_angle=para_in{ni}.rot_angle;
            for i=angle_list
                option1.size_list=size_list;
                option1.comp_size=comp_size_thickness;
                option1.crop_size=para_in{ni}.crop_size;
                option1.rotation_angle_out=i;
                [ img_align1 ] = align_PACBED( database{ni}, [], option1 );
                for j=1:length(size_list)
                    if mask_width>50
                        intensity_mask=0;
                    else
                        img_test0_inv=double(img_align1{j,1}(:,:,1)).*rm_matrix_inv;
                        img_test0_inv(isnan(img_test0_inv))=[];
                        intensity_mask=mean(img_test0_inv);
                    end
                    img_test0=double(img_align1{j,1}(:,:,1)).*rm_matrix;
                    img_test0(isnan(img_test0))=intensity_mask;
                    img_align1{j,1}=cat(3,img_test0,img_test0,img_test0);
                    num=num+1;
                    big_database(:,:,:,num)=img_align1{j,1};
                    if i==0 && size_list(j)==0
                        database_tilt_align(:,:,:,ni)=img_align1{j,1};
                    end
                end
            end
        end
        big_database_tilt=big_database;
    else
        for ni=1:total_img_count
            option1.center=para_in{ni}.center;
            option1.rotation_angle=para_in{ni}.rot_angle;
            for i=angle_list
                option1.size_list=size_list;
                option1.comp_size=[comp_size_thickness,comp_size_tilt];
                option1.crop_size=para_in{ni}.crop_size;
                option1.rotation_angle_out=i;
                [ img_align1 ] = align_PACBED( database{ni}, [], option1 );
                for j=1:length(size_list)
                    for k=1:2
                        if mask_width>50
                            intensity_mask=0;
                        else
                            img_test0_inv=double(img_align1{j,k}(:,:,1)).*rm_matrix_inv;
                            img_test0_inv(isnan(img_test0_inv))=[];
                            intensity_mask=mean(img_test0_inv);
                        end
                        img_test0=double(img_align1{j,k}(:,:,1)).*rm_matrix;
                        img_test0(isnan(img_test0))=intensity_mask;
                        img_align1{j,k}=cat(3,img_test0,img_test0,img_test0);
                    end
                    num=num+1;
                    big_database(:,:,:,num)=img_align1{j,1};
                    big_database_tilt(:,:,:,num)=img_align1{j,2};
                    if i==0 && size_list(j)==0
                        database_tilt_align(:,:,:,ni)=img_align1{j,2};
                    end
                end
            end
        end
    end

    if isempty(angle_list(angle_list==0)) || isempty(size_list(size_list==0)) %input not contain 0 degree/tilt
        for ni=1:total_img_count
            option11.mode='process_only_single';
            option11.center=para_in{ni}.center;
            option11.rotation_angle=para_in{ni}.rot_angle;
            option11.crop_size=para_in{ni}.crop_size*(1+comp_size_thickness/100);
            option11.rotation_angle_out=0;
            [ img_align2 ] = align_PACBED( database{ni}, [], option11 );
            if mask_width>50
                intensity_mask=0;
            else
                img_test0_inv=double(img_align2(:,:,1)).*rm_matrix_inv;
                img_test0_inv(isnan(img_test0_inv))=[];
                intensity_mask=mean(img_test0_inv);
            end
            img_test0=double(img_align2(:,:,1)).*rm_matrix;
            img_test0(isnan(img_test0))=intensity_mask;
            img_align2=cat(3,img_test0,img_test0,img_test0);
            database_tilt_align(:,:,:,ni)=img_align2;
        end
    end    
end
%%data compensation for tilt measurement
total_count=size(big_database,4);
for i=1:total_count
    img_test0=big_database(:,:,:,i);
    big_database(:,:,:,i)=img_test0-min(min(min(img_test0)));
end

%%table list prepare
Class_table_thickness=net_str.convnet_thickness.Layers(end,1).ClassNames;
Score_table_thickness=zeros(length(Class_table_thickness),1);
for i=1:length(Class_table_thickness)
    Score_table_thickness(i,1)=str2double(Class_table_thickness{i,1}(5:end-2));
end

Class_table_tilt=net_str.convnet_tilt.Layers(end,1).ClassNames;
tilt_angle_table=cell(length(Class_table_tilt),1);
tilt_r_table=cell(length(Class_table_tilt),1);
for i=1:length(Class_table_tilt)
    a=Class_table_tilt{i,1};
    num=regexp(a, '_');
    if length(num)>1
        a=a(num(1)+1:num(2)-1);
    else
        a=a(num+1:end);
    end
    num=regexp(a, 'G');
    str1=a(2:num-1);
    str1=strrep(str1,'D','.');
    str2=a(num+1:end);
    str2=strrep(str2,'D','.');
    tilt_angle_table{i}=[str2double(str1),str2double(str2)];
    tilt_r_table{i}.Spec_tilt=sqrt(tilt_angle_table{i}(1)^2+tilt_angle_table{i}(2)^2); %mrad
    tilt_r_table{i}.Spec_azimuth=atand(tilt_angle_table{i}(2)/tilt_angle_table{i}(1)); %degree
end

%%CNN prediction
pred_data_comp=big_database;
scores_t0=zeros(total_count,1);
for iter=1:3
    if cal_gpu==1
        [~,scores] = classify(net_str.convnet_thickness,pred_data_comp);
    else
        [~,scores] = classify(net_str.convnet_thickness,pred_data_comp, 'ExecutionEnvironment','cpu');
    end
    
    scores_t=scores*Score_table_thickness;
    Int_comp=uint8(scores_t*thickness_comp_ratio);
    
    scores_t_diff=abs(scores_t-scores_t0);
    if max(scores_t_diff)<3
        break;
    end
    
    for i=1:total_count
        pred_data_comp(:,:,:,i)=big_database(:,:,:,i)-Int_comp(i);
    end
    scores_t0=scores_t;
end

if cal_gpu==1
    [~,scores_tilt] = classify(net_str.convnet_tilt,big_database_tilt);
else
    [~,scores_tilt] = classify(net_str.convnet_tilt,big_database_tilt, 'ExecutionEnvironment','cpu');
end

scores_red=zeros(ni,size(scores,2));
scores_red_tilt=zeros(ni,size(scores_tilt,2));
block_size=total_count/total_img_count;
for i=1:total_img_count
    scores_red(i,:)=sum(scores((i-1)*block_size+1:i*block_size,:),1)/block_size;
    scores_red_tilt(i,:)=sum(scores_tilt((i-1)*block_size+1:i*block_size,:),1)/block_size;
end
[~,pred_t]=max(scores_red,[],2);
pred_t(:)=Score_table_thickness(pred_t(:));
pred_tilt=cell(ni,1);
pred_tilt_r=zeros(ni,1);
[~,pred_tilt_M]=max(scores_red_tilt,[],2);
for i=1:ni
    pred_tilt{i}=tilt_angle_table{pred_tilt_M(i)};
    pred_tilt_r(i)=tilt_r_table{pred_tilt_M(i)}.Spec_tilt;
end

%%Data analysis
prob_thickness=zeros(total_img_count,1);
prob_tilt=zeros(total_img_count,1);
data_out_all=cell(total_img_count,1);
for i=1:total_img_count 
    scores_pred=scores_red(i,:);
    scores_pred=[scores_pred; 1:size(scores_red,2)]';
    scores_pred=sortrows(scores_pred,-1);
    scores_pred_tilt=scores_red_tilt(i,:);
    scores_pred_tilt=[scores_pred_tilt; 1:size(scores_red_tilt,2)]';
    scores_pred_tilt=sortrows(scores_pred_tilt,-1);
    prob_thickness(i,1)=scores_pred(1,1)*100;
    prob_tilt(i,1)=scores_pred_tilt(1,1)*100;
    data_out_all{i,1}.scores_pred=scores_pred;
    data_out_all{i,1}.scores=scores;
    data_out_all{i,1}.Score_table_thickness=Score_table_thickness;
    data_out_all{i,1}.scores_pred_tilt=scores_pred_tilt;
    data_out_all{i,1}.scores_tilt=scores_tilt;
    data_out_all{i,1}.tilt_angle_table=tilt_angle_table;
    data_out_all{i,1}.tilt_r_table=tilt_r_table;
end

if symmetry_rot_offset~=0
    for i=1:size(database_tilt_align,4)
        database_tilt_align(:,:,:,i)=imrotate(database_tilt_align(:,:,:,i),symmetry_rot_offset,'crop');
    end
end

[ tilt_HG, tilt_azimuth, tilt_align_cell, chk_mark ] = tilt_convert_full_coor_inline( database_tilt_align, pred_tilt, tilt_symmetry );

if chk_reshape==1
    if thickness_stable~=0
        t_mid=median(pred_t);
        for i=1:length(pred_t)
            if abs(pred_t(i)-t_mid)>thickness_stable %out of range
               scores_pred=data_out_all{i,1}.scores_pred;
               for search_num=2:length(scores_pred)
                   t=Score_table_thickness(scores_pred(search_num,2));
                   if abs(t-t_mid)<=thickness_stable
                       pred_t(i)=t;
                       prob_thickness(i)=scores_pred(search_num,1)*100;
                       break;
                   end
               end
            end
        end
    end
    thickness_determine=reshape(pred_t,[size_input(2) size_input(1)])';
    prob_thickness=reshape(prob_thickness,[size_input(2) size_input(1)])';
    thickness_determine=cat(3,thickness_determine,prob_thickness);
    tilt_H=reshape(tilt_HG(:,1),[size_input(2) size_input(1)])';
    tilt_G=reshape(tilt_HG(:,2),[size_input(2) size_input(1)])';
    tilt_azimuth=reshape(tilt_azimuth,[size_input(2) size_input(1)])';
    pred_tilt_r=reshape(pred_tilt_r,[size_input(2) size_input(1)])';
    chk_mark=reshape(chk_mark,[size_input(2) size_input(1)])';
    prob_tilt=reshape(prob_tilt,[size_input(2) size_input(1)])';
    tilt_determine=cat(3,tilt_H,tilt_G,pred_tilt_r,tilt_azimuth,prob_tilt,chk_mark);   
    data_out_all=reshape(data_out_all,[size_input(2) size_input(1)])';
    tilt_align_cell=reshape(tilt_align_cell,[size_input(2) size_input(1)])';
else
    thickness_determine=[pred_t prob_thickness];
    tilt_determine=[tilt_HG,pred_tilt_r,tilt_azimuth,prob_tilt,chk_mark];
end

end



function [ tilt_HG, tilt_azimuth, img_align_tilt, chk_mark ] = tilt_convert_full_coor_inline( img_align, pred_tilt, tilt_symmetry )
%Get full tilt info in four quardrant coor and in spherical coordinate
%Weizong Xu, May, 2017

%analysis four quadrant intensity
tilt_HG=zeros(length(pred_tilt),2);
for i=1:length(pred_tilt)
    tilt_HG(i,1)=pred_tilt{i}(1);
    tilt_HG(i,2)=pred_tilt{i}(2);
end

chk_mark=zeros(size(img_align,4),1);
img_align_tilt=cell(size(img_align,4),1);
ismuSTEM=0;

size_img=size(img_align(:,:,1,1));
center_a=floor((size_img+1)/2);

%generate mask
% mask_III=ones(center_a(1)-1,center_a(2)-1,'uint8');
mask_III=tril(ones(center_a(1)-1,center_a(2)-1,'uint8'),round(center_a(1)/3));
mask_IV=rot90(mask_III);
mask_I=rot90(mask_IV);
mask_II=rot90(mask_I);

for i=1:size(img_align,4)
%     disp(['#',num2str(i)])
    a=img_align(:,:,1,i);
    a_I=a(1:center_a(1)-1,center_a(2)+1:end);
    a_II=a(1:center_a(1)-1,1:center_a(2)-1);
    a_III=a(center_a(1)+1:end,1:center_a(2)-1);
    a_IV=a(center_a(1)+1:end,center_a(2)+1:end);
    int_I=sum(sum(a_I.*mask_I));
    int_II=sum(sum(a_II.*mask_II));
    int_III=sum(sum(a_III.*mask_III));
    int_IV=sum(sum(a_IV.*mask_IV));
    
    if tilt_symmetry==4
        int_quardrant(1,1)=int_I;
        int_quardrant(1,2)=1;
        int_quardrant(2,1)=int_II;
        int_quardrant(2,2)=2;
        int_quardrant(3,1)=int_III;
        int_quardrant(3,2)=3;
        int_quardrant(4,1)=int_IV;
        int_quardrant(4,2)=4;
        int_quardrant=sortrows(int_quardrant,-1);
    
        if int_quardrant(1,2)==1 && int_quardrant(2,2)==2
            a=flip(a,2); %flip from left to right
            tilt_HG(i,1)=-tilt_HG(i,1);
        end

        if int_quardrant(1,2)==1 && int_quardrant(2,2)==3
    %         disp('Very strange pattern! Please check')
            chk_mark(i)=1;
            a=rot90(a);
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=-tilt_HG(i,2);
            tilt_HG(i,2)=tmpH;
        end

        if int_quardrant(1,2)==1 && int_quardrant(2,2)==4
            a=rot90(a);
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=-tilt_HG(i,2);
            tilt_HG(i,2)=tmpH;
        end

        if int_quardrant(1,2)==2 && int_quardrant(2,2)==3
            a=a'; %inverse a
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=tilt_HG(i,2);
            tilt_HG(i,2)=tmpH;
        end

        if int_quardrant(1,2)==2 && int_quardrant(2,2)==4
    %         disp('Very strange pattern! Please check')
            chk_mark(i)=1;
        end

        if int_quardrant(1,2)==3 && int_quardrant(2,2)==1
    %         disp('Very strange pattern! Please check')
            chk_mark(i)=1;
            a=rot90(a,-1);
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=tilt_HG(i,2);
            tilt_HG(i,2)=-tmpH;
        end

        if int_quardrant(1,2)==3 && int_quardrant(2,2)==2
            a=rot90(a,-1);
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=tilt_HG(i,2);
            tilt_HG(i,2)=-tmpH;
        end

        if int_quardrant(1,2)==3 && int_quardrant(2,2)==4
            a=rot90(a,-1)';
            tilt_HG(i,2)=-tilt_HG(i,2);
        end

        if int_quardrant(1,2)==4 && int_quardrant(2,2)==1
            a=rot90(a,2)';
            tmpH= tilt_HG(i,1);
            tilt_HG(i,1)=-tilt_HG(i,2);
            tilt_HG(i,2)=-tmpH;
        end

        if int_quardrant(1,2)==4 && int_quardrant(2,2)==2
    %         disp('Very strange pattern! Please check')
            chk_mark(i)=1;
            a=rot90(a,2);
            tilt_HG(i,1)=-tilt_HG(i,1);
            tilt_HG(i,2)=-tilt_HG(i,2);
        end

        if int_quardrant(1,2)==4 && int_quardrant(2,2)==3
            a=rot90(a,2);
            tilt_HG(i,1)=-tilt_HG(i,1);
            tilt_HG(i,2)=-tilt_HG(i,2);
        end
    end

    if tilt_symmetry==2
        int_quardrant_12=int_I+int_II;
        int_quardrant_23=int_II+int_III;
        int_quardrant_34=int_III+int_IV;
        int_quardrant_41=int_IV+int_I;
        
        if (int_quardrant_12<int_quardrant_34) && (int_quardrant_23<int_quardrant_41)
            a=rot90(a,2);
            tilt_HG(i,1)=-tilt_HG(i,1);
            tilt_HG(i,2)=-tilt_HG(i,2);
        end
        
        if (int_quardrant_12<int_quardrant_34) && (int_quardrant_23>=int_quardrant_41)
            a=flip(a,1); %flip from up to down
            tilt_HG(i,2)=-tilt_HG(i,2);
        end
        
        if (int_quardrant_12>=int_quardrant_34) && (int_quardrant_23<int_quardrant_41)
            a=flip(a,2); %flip from left to right
            tilt_HG(i,1)=-tilt_HG(i,1);
        end
    end    
    
    if ismuSTEM==1 %180 deg rotation of std data of muSTEM relative to mbfit
        a=rot90(a,2);
    end
    
    img_align_tilt{i}=a;

end

%cal zaimuth angle
tilt_azimuth=zeros(size(tilt_HG,1),1);
for i=1:size(tilt_HG,1)
    if tilt_HG(i,2)>0 && tilt_HG(i,1)<0 %I
        tilt_azimuth(i)=atand(-tilt_HG(i,2)/tilt_HG(i,1)); %degree
    end
    
    if tilt_HG(i,2)>0 && tilt_HG(i,1)>0 %II
        tilt_azimuth(i)=180-atand(tilt_HG(i,2)/tilt_HG(i,1)); %degree
    end
    
    if tilt_HG(i,2)<0 && tilt_HG(i,1)>0 %III
        tilt_azimuth(i)=atand(-tilt_HG(i,2)/tilt_HG(i,1))+180; %degree
    end
    
    if tilt_HG(i,2)<0 && tilt_HG(i,1)<0 %IV
        tilt_azimuth(i)=360-atand(tilt_HG(i,2)/tilt_HG(i,1)); %degree
    end

    if tilt_HG(i,2)==0 && tilt_HG(i,1)<=0 %
        tilt_azimuth(i)=0; %degree
    end
    
    if tilt_HG(i,2)==0 && tilt_HG(i,1)>0
        tilt_azimuth(i)=180; %degree
    end
    
    if tilt_HG(i,2)>0 && tilt_HG(i,1)==0
        tilt_azimuth(i)=90; %degree
    end
    
    if tilt_HG(i,2)<0 && tilt_HG(i,1)==0
        tilt_azimuth(i)=270; %degree
    end
end

end