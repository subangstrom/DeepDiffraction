function Img_series_Looper(imageSeries,name_list)
%2D version of series loop, convenient to check image database
%Weizong Xu, April, 2017

if ~exist('name_list','var')
    name_list=[];
end

if isempty(imageSeries)
    disp('Empty input.')
    return;
end

chk1(1)=size(imageSeries,1);
chk1(2)=size(imageSeries,2);
if chk1(1)<chk1(2)
    imageSeries=imageSeries';
    chk1=flip(chk1,2);
end

chk2(1)=size(name_list,1);
chk2(2)=size(name_list,2);
if chk2(1)<chk2(2)
    name_list=name_list';
    chk2=flip(chk2,2);
end

if abs(sum(chk1-chk2))>0 && ~isempty(name_list)
    disp('Input database and its description are not consistent. Only show image from database!')
    name_list=[];
end

global currentFrame_x  currentFrame_y
currentFrame_x = 1;
currentFrame_y = 1;
currentFrame_s = 1;

for i=1:length(imageSeries)
    imageSeries{i}=imageSeries{i}(:,:,1); %rgb to gray if have
end

if isa(imageSeries{1},'uint8')
    numSteps_min = 256;
    numSteps_max = 256;
else
    if isa(imageSeries{1},'uint16')
        numSteps_min = 65536;
        numSteps_max = 65536;
    else
        for i=1:size(imageSeries,1)
            for ii=1:size(imageSeries,2)
            imageSeries{i,ii}=uint16((imageSeries{i,ii}-min(min(imageSeries{i,ii})))/(max(max(imageSeries{i,ii}))-min(min(imageSeries{i,ii})))*65536);
            end
        end
            numSteps_min = 65536;
            numSteps_max = 65536;
    end
end

numSteps_x = size(imageSeries,1);
numSteps_y = size(imageSeries,2);
numSteps_s = round(min(size(imageSeries{1},1),size(imageSeries{1},1))/2);
Int_max=max(max(imageSeries{1,1}));
Int_min=min(min(imageSeries{1,1}));

figure('color',[1,1,1]*0.92)


imagesc(imageSeries{currentFrame_x,currentFrame_y},[Int_min, Int_max]);
if ~isempty(name_list)
    title(strrep(name_list{currentFrame_x,currentFrame_y},'_','-'))
end
colormap(gray);
axis image;

if numSteps_x>1
    h_x = uicontrol('Style', 'slider', 'Min',1,'Max',numSteps_x,'Value',currentFrame_x,'Callback', @setImage);
    set(h_x, 'SliderStep', [1/(numSteps_x-1) , 1/(numSteps_x-1) ]);
    set(h_x,  'Position', [135 0 300 25]);
end

if numSteps_y>1
    h_y = uicontrol('Style', 'slider', 'Min',1,'Max',numSteps_y,'Value',currentFrame_y,'Callback', @setImage);
    set(h_y, 'SliderStep', [1/(numSteps_y-1) , 1/(numSteps_y-1) ]);
    set(h_y,  'Position', [50 70 25 300]);
end

h_s = uicontrol('Style', 'slider', 'Min',1,'Max',numSteps_s,'Value',currentFrame_s,'Callback', @setImage);
set(h_s, 'SliderStep', [1/(numSteps_s-1) , 1/(numSteps_s-1) ]);
set(h_s,  'Position', [470 70 25 300]);

h_min = uicontrol('Style', 'slider', 'Min',-numSteps_max+1,'Max',numSteps_max-1,'Value',Int_min,'Callback', @setImage);
set(h_min, 'SliderStep', [1/(numSteps_min-1)/2 , 1/(numSteps_min-1)/2 ]);
set(h_min,  'Position', [500 70 25 300]);

h_max = uicontrol('Style', 'slider', 'Min',-numSteps_max+1,'Max',numSteps_max-1,'Value',Int_max,'Callback', @setImage);
set(h_max, 'SliderStep', [1/(numSteps_max-1)/2 , 1/(numSteps_max-1)/2 ]);
set(h_max,  'Position', [530 70 25 300]);

if numSteps_x>1 && numSteps_y>1
    text = uicontrol('Style','text','String','Frame: 1 1');
else
    text = uicontrol('Style','text','String','Frame: 1');
end
set(text, 'Position', [470 370 80 20]);

text_size = uicontrol('Style','text','String','x1');
set(text_size, 'Position', [470 50 30 20]);

text_min = uicontrol('Style','text','String',num2str(Int_min));
set(text_min, 'Position', [500 50 30 20]);

text_max = uicontrol('Style','text','String',num2str(Int_max));
set(text_max, 'Position', [530 50 30 20]);

str_x='Loop_x';
str_y='Loop_y';
if numSteps_x==1 || numSteps_y==1
    str_x='Loop';
    str_y='Loop';
end

if numSteps_x>1
    btn_x = uicontrol('Style', 'pushbutton', 'String', str_x,...
        'Position', [455 5 50 20], 'Interruptible','on', ...
        'Callback', @loop_x);
end

if numSteps_y>1
    btn_y = uicontrol('Style', 'pushbutton', 'String', str_y,...
        'Position', [505 5 50 20], 'Interruptible','on', ...
        'Callback', @loop_y);
end

btn_reset = uicontrol('Style', 'pushbutton', 'String', 'Reset',...
    'Position', [490 30 50 20], 'Interruptible','on', ...
    'Callback', @reset_int);

isPlaying = false; %func when press stop
x_run=0;
y_run=0;

    function setImage(event_obj, var)
        if numSteps_x==1
            currentFrame_x=1;
        else
            currentFrame_x = round(get(h_x,'value'));
        end
        
        if numSteps_y==1
            currentFrame_y=1;
        else
            currentFrame_y = round(get(h_y,'value'));
        end
        currentFrame_s = round(get(h_s,'value'));
        Int_min = round(get(h_min,'value'));
        Int_max = round(get(h_max,'value'));
        if currentFrame_x>numSteps_x
            currentFrame_x=numSteps_x;
        end
        if currentFrame_y>numSteps_y
            currentFrame_y=numSteps_y;
        end
        img_read=imageSeries{currentFrame_x,currentFrame_y};
        img_read=img_read(currentFrame_s:end-currentFrame_s,currentFrame_s:end-currentFrame_s,:);
        if Int_min<Int_max
            imagesc(img_read,[Int_min, Int_max]);
        else
            imagesc(img_read,[Int_max-1, Int_min]);
        end
        if ~isempty(name_list)
            title(strrep(name_list{currentFrame_x,currentFrame_y},'_','-'))
        end
        if numSteps_x==1
            set(text, 'String', ['Frame:', num2str(currentFrame_y)]);
        else
            if numSteps_y==1
                set(text, 'String', ['Frame:', num2str(currentFrame_x)]);
            else
                set(text, 'String', ['Frame:', num2str(currentFrame_x),' ',num2str(currentFrame_y)]);
            end
        end
        set(text_size, 'String', ['x', num2str(round(numSteps_s/(numSteps_s-currentFrame_s)*100)/100)]);
        set(text_min, 'String', num2str(round(Int_min)));
        set(text_max, 'String', num2str(round(Int_max)));
        colormap(gray);
        axis image;
    end


    function loop_x(event_obj, var)
        
        if isPlaying == false || y_run==1
            isPlaying = true;
            set(btn_x, 'String', 'stop');
            if numSteps_y > 1
                set(btn_y, 'String', str_y);
            end
            x_run=1;
            y_run=0;
        else
            isPlaying = false;
            set(btn_x, 'String', str_x);
            x_run=0;
            y_run=0;
        end
        
        while(isPlaying)
            currentFrame_x = currentFrame_x + 1;
            set(h_x,  'Value', currentFrame_x);
            pause(1/20);
            setImage(event_obj, var)
            if currentFrame_x == numSteps_x
                currentFrame_x = 1;
            end
        end
    end

    function loop_y(event_obj, var)
        
        if isPlaying == false || x_run==1
            isPlaying = true;
            set(btn_y, 'String', 'stop');
            if numSteps_x > 1
                set(btn_x, 'String', str_x);
            end
            x_run=0;
            y_run=1;
        else
            isPlaying = false;
            set(btn_y, 'String', str_y);
            x_run=0;
            y_run=0;
        end        

        while(isPlaying)
            currentFrame_y = currentFrame_y + 1;
            set(h_y,  'Value', currentFrame_y);
            pause(1/20);
            setImage(event_obj, var)
            if currentFrame_y == numSteps_y
                currentFrame_y = 1;
            end
        end
    end

function reset_int(event_obj, var)
    Int_max=max(max(imageSeries{1,1}));
    Int_min=min(min(imageSeries{1,1}));
    currentFrame_s=1;
    set(h_s,'value',1);
    set(h_min,'value',Int_min);
    set(h_max,'value',Int_max);
    setImage(event_obj, var)
end
end

