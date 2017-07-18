function [ tilt_HG, tilt_r, tilt_azimuth, img_align_cell_tilt, chk_mark ] = tilt_convert_full_coor( img_align_cell, tilt_determine, tilt_symmetry )
%Get full tilt info in four quardrant coor and in spherical coordinate
%Weizong Xu, April, 2017

%datainput must be 2D
%analysis four quadrant intensity

if ~exist('tilt_symmetry','var')
    tilt_symmetry=4;
    disp('Set symmetry as four-fold.')
end
if tilt_symmetry~=4 && tilt_symmetry~=2
    disp('The symmetry setting is currently not supported. Auto change it to 2-fold symmetry to continue the processing.')
    tilt_symmetry=2;
end

tilt_r=zeros(size(tilt_determine,1),1);
tilt_HG=zeros(size(tilt_determine,1),2);
for i=1:size(tilt_determine,1)
    tilt_r(i,1)=norm(tilt_determine{i,1});
    tilt_HG(i,1)=tilt_determine{i,1}(1);
    tilt_HG(i,2)=tilt_determine{i,1}(2);
end

chk_mark=zeros(length(img_align_cell),1);
ismuSTEM=0;

size_img=size(img_align_cell{1}(:,:,1));
center_a=floor((size_img+1)/2);

%generate mask
% mask_III=ones(center_a(1)-1,center_a(2)-1,'uint8');
mask_III=tril(ones(center_a(1)-1,center_a(2)-1,'uint8'),round(center_a(1)/3));
mask_IV=rot90(mask_III);
mask_I=rot90(mask_IV);
mask_II=rot90(mask_I);

for i=1:length(img_align_cell)
    disp(['#',num2str(i)])
    a=img_align_cell{i,1};
    a=a(:,:,1);
    a_I=a(1:center_a(1)-1,center_a(2)+1:end);
    a_II=a(1:center_a(1)-1,1:center_a(2)-1);
    a_III=a(center_a(1)+1:end,1:center_a(2)-1);
    a_IV=a(center_a(1)+1:end,center_a(2)+1:end);
    int_I=sum(sum(a_I.*mask_I));
    int_II=sum(sum(a_II.*mask_II));
    int_III=sum(sum(a_III.*mask_III));
    int_IV=sum(sum(a_IV.*mask_IV));
    
% %     int_quardrant(1,1)=sum(sum(a_I.*mask_I));
% %     int_quardrant(1,2)=1;
% %     int_quardrant(2,1)=sum(sum(a_II.*mask_II));
% %     int_quardrant(2,2)=2;
% %     int_quardrant(3,1)=sum(sum(a_III.*mask_III));
% %     int_quardrant(3,2)=3;
% %     int_quardrant(4,1)=sum(sum(a_IV.*mask_IV));
% %     int_quardrant(4,2)=4;
% %     int_quardrant=sortrows(int_quardrant,-1);
    
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
            disp('Very strange pattern! Please check')
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
            disp('Very strange pattern! Please check')
            chk_mark(i)=1;
        end

        if int_quardrant(1,2)==3 && int_quardrant(2,2)==1
            disp('Very strange pattern! Please check')
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
            disp('Very strange pattern! Please check')
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
    % figure;imshow(a,[])
    if ismuSTEM==1 %180 deg rotation of std data of muSTEM relative to mbfit
        a=rot90(a,2);
    end
    
    img_align_cell_tilt{i,1}=a;
end
% Img_series_Looper(img_align_cell_tilt,calb_file)
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

