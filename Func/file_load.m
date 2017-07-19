function [ database, name_list ] = file_load( filepath, subfolder )
%load file to database for analysis
%Weizong Xu, July, 2017

database={};
name_list={};

if ~exist('filepath','var') || isempty(filepath)
	[name_list, pathname] = uigetfile( ...
     {'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.tiff', 'All image Files (*.jpg, *.jpeg, *.png, *.bmp, *.tif, *.tiff)';
      '*.jpg;*.jpeg',  'JPEG format (*.jpg, *.jpeg)'; ...
      '*.png','PNG format (*.png)'; ...
      '*.bmp','BMP format (*.bmp)'; ...
      '*.tif;*.tiff','TIFF format (*.tif, *.tiff)'; ...
      '*.*',  'All Files (*.*)'}, ...
      'Select a file or a series of files','MultiSelect','on');

    if ~iscell(name_list)
        if name_list==0
            error('No file is selected.')
        end
        name_list={name_list};
    end
    name_list=name_list';
    database=cell(length(name_list),1);
    for i=1:length(name_list)
        if ispc
            database{i}=imread([pathname,'\',name_list{i}]);
        else
            database{i}=imread([pathname,'/',name_list{i}]);
        end
        database{i}=database{i}(:,:,1);
    end
else
    if ~exist('subfolder','var') || subfolder~=1
    dirData = dir(filepath);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    for i=1:length(fileList)
        if ispc
            fileList{i}=[filepath,'\',fileList{i}];
        else
            fileList{i}=[filepath,'/',fileList{i}];
        end
    end
    else
        fileList = getAllFiles(filepath);
    end
    if isempty(fileList)
        error('No file found in this folder, please check!')
    end
    num=0;
    for i=1:length(fileList)
        [~,name,~]=fileparts(fileList{i});
        try
            database{num+1}=imread(fileList{i});
            database{num+1}=database{num+1}(:,:,1);
            name_list{num+1}=name;
            num=num+1;
        catch
            %ignore this image
        end
    end
end

if isempty(name_list)
    error('No image file found in this folder, please check!')
end

end

function fileList = getAllFiles(dirName)

  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
  end

end