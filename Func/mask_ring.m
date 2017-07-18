function [ img_in, img_inv ] = mask_ring( img_in, r_range, set_value )
%set zero value on image at a radius range;
%note: this func is slow, generate a mask from this func, not use in for loop
%Weizong Xu, May, 2017

if ~exist('r_range','var')
    return;
end

if ~exist('set_value','var')
    set_value=0;
end

ring_center=[(size(img_in,1)+1)/2 (size(img_in,2)+1)/2];
for i=1:size(img_in,1)
    for j=1:size(img_in,2)
        dist=sqrt((i-ring_center(1))^2+(j-ring_center(2))^2);
        if dist>min(r_range) && dist<max(r_range)
            img_in(i,j)=set_value;
        end
    end
end

%produce inverse of the mask
img_inv=img_in;
img_inv(isnan(img_inv))=0;
img_inv(img_inv==1)=nan;
img_inv(img_inv==0)=1;

end

