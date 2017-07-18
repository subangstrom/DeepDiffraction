function [ A2, B_crop, AB_diff ] = func_lsf_PACBED_evaluate( img_in, PACBED_data, LSF_opt, tilt_symmetry )
%Evaluate LSF matching quality
%Weizong Xu, May, 2017

if ~exist('tilt_symmetry','var')
    tilt_symmetry=4;
    disp('Set symmetry as four-fold.')
end
if tilt_symmetry~=4 && tilt_symmetry~=2
    disp('The symmetry setting is currently not supported. Auto change it to 2-fold symmetry to continue the processing.')
    tilt_symmetry=2;
end

rot_angle=LSF_opt(2);
img_size=LSF_opt(3);
int_comp=LSF_opt(4);
t_select=LSF_opt(5);
i_offsetx=LSF_opt(6);
i_offsety=LSF_opt(7);
tilt_select=LSF_opt(8);
num_iter=LSF_opt(9);

A20=double(img_in);
A20=A20/max(max(A20))*255;%normalize image back to 0-255 scale
A2=A20;

if num_iter==2 && tilt_symmetry==2 %2-fold symmetry full check
    A2=rot90(A20,2);
end
if num_iter==3 && tilt_symmetry==2
    A2=flip(A20,1);
end
if num_iter==4 && tilt_symmetry==2
    A2=flip(A20,2);
end
    
if num_iter==2 && tilt_symmetry==4
    A2=A20';
end
if num_iter==3 && tilt_symmetry==4
    A2=rot90(A20,2);
end
if num_iter==4 && tilt_symmetry==4
    A2=rot90(A20,2)';
end
if num_iter==5
    A2=flip(A20,2);
end
if num_iter==6
    A2=rot90(A20);
end
if num_iter==7
    A2=rot90(A20,-1);
end
if num_iter==8
    A2=rot90(A20,-1)';
end

B_load=PACBED_data{tilt_select,t_select};
B_load=double(B_load(:,:,1));
B_load=B_load/max(max(B_load))*255;%normalize image back to 0-255 scale
B_rot=imrotate(B_load,rot_angle,'crop');
B=imresize(B_rot,img_size);
B_center0=round(size(B)/2);
B_center=[B_center0(1)+i_offsetx, B_center0(2)+i_offsety];
B_crop=B(B_center(1)-113:B_center(1)+113,B_center(2)-113:B_center(2)+113);
B_crop=double(B_crop);
AB_diff=A2-B_crop*int_comp;
% AB_diff=AB_diff.^2;


end

