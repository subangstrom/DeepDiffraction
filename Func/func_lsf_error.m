function [ data_out ] = func_lsf_error( img_in, PACBED_data, s_list, threshold )
%Least Square method to determine PACBED thickness and mistilt
%Weizong Xu, June, 2017

img_size=s_list.img_size;
int_comp=s_list.int_comp;
tilt_select_list=s_list.tilt_select_list;
t_select_list=s_list.t_select_list;
name_list=s_list.name_list;
filename=s_list.filename;

data_out=cell(1,5);
A2=double(img_in);

%%main
num=0;AB_list=zeros(1,9);AB_list2=zeros(1,9);tmp_max=1e15;
    
for i_thickness=1:length(t_select_list)
    t_select=t_select_list(i_thickness);
    for j_tilt=1:length(tilt_select_list)
        tilt_select=tilt_select_list(j_tilt);
%         disp(['Searching ',name_list{tilt_select,t_select},' for ',filename])
        B_load=PACBED_data{tilt_select,t_select};
        B_load=double(B_load(:,:,1));
        B_load=B_load/max(max(B_load))*255;%normalize image back to 0-255 scale
        B=imresize(B_load,img_size);
        B_center=round(size(B)/2);
        B_crop=B(B_center(1)-113:B_center(1)+113,B_center(2)-113:B_center(2)+113);
        B_crop=double(B_crop);
        AB_diff=A2-B_crop*int_comp;
        AB_diff=AB_diff.^2;
        AB_diff_max=sum(sum(AB_diff));
        if AB_diff_max<=threshold && length(AB_list2)<5000
            num=num+1;
            AB_list2(num,:)=[AB_diff_max, 0, img_size, int_comp, t_select, 0, 0, tilt_select, 1];
        end
        if AB_diff_max<tmp_max
            AB_list=[AB_diff_max, 0, img_size, int_comp, t_select, 0, 0, tilt_select, 1];
            tmp_max=AB_diff_max;
            if AB_diff_max<=threshold && length(AB_list2)>=5000
                num=num+1;
                AB_list2(num,:)=AB_list;
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

