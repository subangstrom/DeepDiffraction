function [ img_series, background ] = image_series_BG_rm( img_series, type, klim )
%DFT method to remove background
%klim=1, just simple background removal
%Weizong Xu

if ~exist('type','var')
    %simple remove background around 4 corners
    type=0;
end

if ~exist('klim','var')
    klim=1;
end

if ~iscell(img_series)
    img_series={img_series};
    chk_cell=0;
else
    chk_cell=1;
end

background=cell(size(img_series,1),size(img_series,2));
for i=1:size(img_series,1)
    if size(img_series,1)*size(img_series,2)>1000
        disp(['Processing #',num2str(i)])
    end
    for j=1:size(img_series,2)
        image=double(img_series{i,j});
        image0=image;
        if size(image,3)==3
            chk_size=3;
            image=image(:,:,1);
        else
            chk_size=1;
        end
        if type==1
            [image, background{i,j}] = clearImageBackground(image, klim);
        end
        image=uint8(image);
        image0=uint8(image0);
        image=image-max([image(1,1),image(1,end),image(end,1),image(end,end)])+1;
        if type~=1
            background{i,j}=image0-image;
        end
        image=double(image);
        image=image/max(max(image))*255; %normalize image for NN
        
        if chk_size==3
            image=cat(3,image,image,image);
        end
        img_series{i,j}=uint8(image);
    end
end

if chk_cell==0
    img_series=img_series{1};
end

end

