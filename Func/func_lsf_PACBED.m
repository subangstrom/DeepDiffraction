function [ data_out ] = func_lsf_PACBED( img_in, PACBED_data, s_list, threshold )
%Least Square method to determine PACBED thickness and mistilt
%Weizong Xu, May, 2017

rot_angle_list=s_list.rot_angle_list;
img_size_list=s_list.img_size_list;
int_comp_list=s_list.int_comp_list;
tilt_select_list=s_list.tilt_select_list;
t_select_list=s_list.t_select_list;
name_list=s_list.name_list;
offsetx_list=s_list.offsetx_list;
offsety_list=s_list.offsety_list;
filename=s_list.filename;
chk_mark=s_list.chk_mark;

data_out=cell(1,5);
img_in=img_in-min(min(img_in));
A20=double(img_in);
A20=A20/max(max(A20))*255;%normalize image back to 0-255 scale
if chk_mark==0
    tot_loop=1;
else
    if chk_mark==1 || chk_mark==-1
        tot_loop=4;
    else
        tot_loop=8;
    end
end
%%main
num=0;AB_list=zeros(1,9);AB_list2=zeros(1,9);tmp_max=1e15;
for num_iter=1:tot_loop
    A2=A20;
    
    if num_iter==2 && chk_mark==-1 %2-fold symmetry full check
        A2=rot90(A20,2);
    end
    if num_iter==3 && chk_mark==-1
        A2=flip(A20,1);
    end
    if num_iter==4 && chk_mark==-1
        A2=flip(A20,2);
    end
    
    if num_iter==2 && chk_mark~=-1 %4-fold symmetry check
        A2=A20';
    end
    if num_iter==3 && chk_mark~=-1
        A2=rot90(A20,2);
    end
    if num_iter==4 && chk_mark~=-1
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
    
    for i_thickness=1:length(t_select_list)
        t_select=t_select_list(i_thickness);
        for j_tilt=1:length(tilt_select_list)
            tilt_select=tilt_select_list(j_tilt);
            disp(['Searching ',name_list{tilt_select,t_select},' for ',filename, ', iter ',num2str(num_iter),'/',num2str(tot_loop)])
            B_load=PACBED_data{tilt_select,t_select};
            B_load=double(B_load(:,:,1));
            B_load=B_load/max(max(B_load))*255;%normalize image back to 0-255 scale
            for k_rot_angle=1:length(rot_angle_list)
                rot_angle=rot_angle_list(k_rot_angle);
                B_rot=imrotate(B_load,rot_angle,'crop');
                for l_img_size=1:length(img_size_list)
                    img_size=img_size_list(l_img_size);
                    B=imresize(B_rot,img_size);
                    B_center0=round(size(B)/2);
                    for o_offsetx=1:length(offsetx_list)
                        i_offsetx=offsetx_list(o_offsetx);
                        for o_offsety=1:length(offsety_list)
                            i_offsety=offsety_list(o_offsety);
                            B_center=[B_center0(1)+i_offsetx, B_center0(2)+i_offsety];
                            B_crop=B(B_center(1)-113:B_center(1)+113,B_center(2)-113:B_center(2)+113);
                            B_crop=double(B_crop);
                            for o_int_comp=1:length(int_comp_list)
                                int_comp=int_comp_list(o_int_comp);
                                AB_diff=A2-B_crop*int_comp;
                                AB_diff=AB_diff.^2;
                                AB_diff_max=sum(sum(AB_diff));
                                if AB_diff_max<=threshold && length(AB_list2)<5000
                                    num=num+1;
                                    AB_list2(num,:)=[AB_diff_max, rot_angle, img_size, int_comp, t_select, i_offsetx, i_offsety, tilt_select, num_iter];
                                end
                                if AB_diff_max<tmp_max
                                    AB_list=[AB_diff_max, rot_angle, img_size, int_comp, t_select, i_offsetx, i_offsety, tilt_select, num_iter];
                                    tmp_max=AB_diff_max;
                                    if AB_diff_max<=threshold && length(AB_list2)>=5000
                                        num=num+1;
                                        AB_list2(num,:)=AB_list;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%%sort
AB_list_sort=sortrows(AB_list);
AB_list_sort2=sortrows(AB_list2);

%%output
t_select=AB_list_sort(1,5);
tilt_select=AB_list_sort(1,8);

data_out{1,1}=AB_list_sort;
data_out{1,2}=filename;
data_out{1,3}=name_list{tilt_select,t_select};
data_out{1,4}=AB_list_sort(5);
data_out{1,5}=AB_list_sort2;


end

