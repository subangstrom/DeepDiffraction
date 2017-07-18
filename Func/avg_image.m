function [ img_out ] = avg_image( img, size_filter )
%Average filter of the image, and avoid dark line at boundary
%Fast algorithm using MATLAB internal function
%Weizong Xu, Dec. 2016

comp_size=floor(size_filter/2);
s_img=size(img);
tmp_img=zeros(s_img(1)+2*comp_size,s_img(2)+2*comp_size);

tmp_img(comp_size+1:end-comp_size,comp_size+1:end-comp_size)=img;

for i=1:comp_size
    tmp_img(i,comp_size+1:end-comp_size)=img(i+1,:);
    tmp_img(end-i+1,comp_size+1:end-comp_size)=img(end-i,:);
    tmp_img(comp_size+1:end-comp_size,i)=img(:,i+1);
    tmp_img(comp_size+1:end-comp_size,end-i+1)=img(:,end-i);
    
end

%deal with corner
tmp_A=img(1:comp_size,1:comp_size);
tmp_img(1:comp_size,1:comp_size)=rot90(tmp_A,2);
tmp_A=img(1:comp_size,end-comp_size:end);
tmp_img(1:comp_size,end-comp_size:end)=rot90(tmp_A,2);
tmp_A=img(end-comp_size:end,1:comp_size);
tmp_img(end-comp_size:end,1:comp_size)=rot90(tmp_A,2);
tmp_A=img(end-comp_size:end,end-comp_size:end);
tmp_img(end-comp_size:end,end-comp_size:end)=rot90(tmp_A,2);


tmp_img_out=filter2(fspecial('average',size_filter),tmp_img);
img_out=tmp_img_out(comp_size+1:end-comp_size,comp_size+1:end-comp_size);
end

