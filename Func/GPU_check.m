function [ cal_gpu ] = GPU_check( convnet )
%check if GPU can be used for NN
%Weizong, May, 2017

cal_gpu=1;
if isempty(convnet)
    disp('Network check warning: file not available.')
    return;
end

try
    gpu_info=gpuDevice;
    if gpu_info.AvailableMemory<2e9
        cal_gpu=0;
    end
catch
    cal_gpu=0;
end


if cal_gpu==1
    try
        img_test0=zeros(convnet.Layers(1,1).InputSize,'uint8')+128;
        classify(convnet,img_test0);
    catch ME
        disp(ME)
        disp('WARNING! GPU is not working properly!')
        cal_gpu=0;
    end
end

end

