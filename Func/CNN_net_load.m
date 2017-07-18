function [ convnet, info ] = CNN_net_load( net_name )
%load CNN network from mat file
%Weizong, April, 2017

try
    net_str=load(net_name);
catch
    disp(['Cannot load ',net_name])
    convnet=[];
    info.error=1;
    return;
end

try
    convnet=net_str.convnetTransfer;
    try
        info=net_str.info;
    catch
        info.error=1;
    end
catch
    convnet=net_str.net;
    info.error=1;
end


end

