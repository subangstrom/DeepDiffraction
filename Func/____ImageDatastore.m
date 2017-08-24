classdef (Sealed) ImageDatastore < ...
                  matlab.io.datastore.CustomReadDatastore & ...
                  matlab.io.datastore.internal.ScalarBase & ...
                  matlab.mixin.CustomDisplay
%The code has been modified for specifit task in the neural network training for Matlab 2017a by Weizong Xu.
%Copy to Matlab installation folder, rename to ImageDatastore.m
%e.g. .../MATLAB_R2017a/toolbox/matlab/datastoreio/+matlab/+io/+datastore/@ImageDatastore/ImageDatastore.m
%For shift and rotation CNN training, comment on the random rotate code portion (line 399-410)

%IMAGEDATASTORE Datastore for a collection of image files.
%   IMDS = imageDatastore(LOCATION) creates an ImageDatastore IMDS given the
%   LOCATION of the image data. LOCATION has the following properties:
%      - Can be a filename or a folder name
%      - Can be a cell array of multiple file or folder names
%      - Can contain a relative path (HDFS requires a full path)
%      - Can contain a wildcard (*) character.
%      - All the files in LOCATION must have the same extension and be
%        supported by IMFORMATS
%
%   IMDS = imageDatastore(__,'IncludeSubfolders',TF) specifies the logical
%   true or false to indicate whether the files in each folder and its
%   subfolders are included recursively or not.
%
%   IMDS = imageDatastore(__,'FileExtensions',EXTENSIONS) specifies the
%   extensions of files to be included. The extensions are not required to
%   be supported by IMFORMATS. Values for EXTENSIONS can be:
%      - A character vector, such as '.jpg' or '.png' (empty quotes '' are
%        allowed for files without extensions)
%      - A cell array of character vectors, such as {'.jpg', '.png'}
%
%   IMDS = imageDatastore(__,'ReadSize',READSIZE) specifies the maximum
%   number of image files to read in a call to the read function. By default,
%   READSIZE is 1. The output of read is a cell array of image data when
%   READSIZE > 1.
%
%   IMDS = imageDatastore(__,'ReadFcn',@MYCUSTOMREADER) specifies the user-
%   defined function to read files. The value of 'ReadFcn' must be a
%   function handle with a signature similar to the following:
%      function data = MYCUSTOMREADER(filename)
%      ..
%      end
%
%   IMDS = imageDatastore(__,'LabelSource',SOURCE) specifies the source from
%   which the Labels property obtains labels. By default, the value of
%   SOURCE is 'none'. If SOURCE is 'foldernames', then the values for the
%   Labels property are obtained from the folder names of the image files.
%
%   IMDS = imageDatastore(__,'Labels',LABELS) specifies the datastore labels
%   according to LABELS. LABELS must be a cell array of character vectors or
%   a vector of numeric, logical, or categorical type.
%
%   ImageDatastore Properties:
%
%      Files           - Cell array of file names
%      ReadSize        - Upper limit on the number of images returned by the read method
%      ReadFcn         - Function handle used to read files
%      Labels          - A set of labels for images
%
%   ImageDatastore Methods:
%
%      hasdata        - Returns true if there is more data in the datastore
%      read           - Reads the next consecutive file
%      reset          - Resets the datastore to the start of the data
%      preview        - Reads the first image from the datastore
%      readimage      - Reads a specified image from the datastore
%      readall        - Reads all image files from the datastore
%      partition      - Returns a new datastore that represents a single
%                       partitioned portion of the original datastore
%      numpartitions  - Returns an estimate for a reasonable number of
%                       partitions to use with the partition function,
%                       according to the total data size
%      splitEachLabel - Splits the ImageDatastore labels according to the
%                       specified proportions, which can be represented as
%                       percentages or number of files.
%      countEachLabel - Counts the number of unique labels in the ImageDatastore
%      shuffle        - Shuffles the files of ImageDatastore using randperm
%
%   Example:
%   --------
%      folders = fullfile(matlabroot,'toolbox','matlab',{'demos','imagesci'});
%      exts = {'.jpg','.png','.tif'};
%      imds = imageDatastore(folders,'FileExtensions',exts);
%      img1 = read(imds);                  % Read the first image
%      img2 = read(imds);                  % Read the next image
%      readall(imds)                       % Read all of the images
%      imgarr = cell(numel(imds.Files),1);
%      for i = 1:numel(imds.Files)         % Read images using a for loop
%          imgarr{i} = readimage(imds,i);
%      end
%
%   See also datastore, mapreduce, imformats, imageDatastore.

%   Copyright 2015-2016 The MathWorks, Inc.

    properties (Dependent)
        %Files -
        % A cell array of image files.
        Files;
    end

    properties (Hidden)
        %CachedRead -
        % 'on' | 'off' - Toggle to use cached read.
        CachedRead = 'on';
        %MaxThreads -
        % numeric - Max number of threads for cached read.
        MaxThreads = 'default';
        % PrefetchSize indicates how much data to read in advance
        PrefetchSize = 0;
    end

    properties
        %ReadSize -
        % Number of image files for one read.
        ReadSize;
        %Labels -
        % A set of labels for images with a one-to-one mapping with Files.
        Labels = {};
    end

    properties (Access = private)
        % deployement needs a way to get files before resolving them
        UnResolvedFiles;
        % On subsasgn'ing Files with [], we need a way to subsasgn Labels
        % to [].
        EmptyIndexes;
        % On subsasgn'ing Files with more files, we need a way to subsasgn
        % Labels to default values.
        AddedIndexes;
        % On assigning Files property with new files, we need to assign
        % default values for Labels property.
        FilesAssigned = false;
        % Whenever we set files, determine number of files instead of
        % calling numel all the time (reset method). We also need to check
        % for changed indexes of Files, to subsasgn Labels.
        % Refer method resolveLabelsIndexes.
        NumFiles;
        % To help support future forward compatibility.  The value
        % indicates the version of MATLAB.
        SchemaVersion;
        % To know if ReadFcn is changed
        IsReadFcnDefault;
        % To know if MaxThreads is changed
        IsMaxThreadsDefault = true;
        % To know if cached read is already used
        IsUsingCachedRead = true;
        % Batch reader
        BatchReader;
        % File extensions
        FileExts;
        % Data buffer holding on to cached images
        DataBuffer;
        % Error buffer pointing to cached images
        ErrorBuffer;
    end

    properties (Constant, Access = private)
        WHOLE_FILE_CUSTOM_READ_SPLITTER_NAME = 'matlab.io.datastore.splitter.WholeFileCustomReadSplitter';
        CONVENIENCE_CONSTRUCTOR_FCN_NAME = 'imageDatastore';
    end
    % Constructor
    methods
        % ImageDataStore can be constructed with files argument, optionally
        % with ReadFcn, IncludeSubfolders, FileExtensions Name-Value pairs.
        function imds = ImageDatastore(files, varargin)
            try
                import matlab.io.datastore.ImageDatastore;
                nv = iParseNameValues(varargin{:});
                initDatastore(imds, files, nv);
                imds.BatchReader = [];
                imds.PrefetchSize = 0;
                imds.UnResolvedFiles = files;
                % SchemaVersion indicates the release number of MATLAB. This will be empty in
                % 14b or the appropriate release, if we set it in the constructor.
                imds.SchemaVersion = version('-release');
            catch e
                throwAsCaller(e);
            end
        end
    end

    % Set and Get methods for properties
    methods
        % Set method for Files
        function set.Files(imds, files)
            try
                c = onCleanup(@() initializeSubsAsgnIndexes(imds));
                [diffIndexes,currIndexes] = setNewFilesAndFileSizes(imds, files);
                resolveLabelsIndexes(imds);
                setFileExtensions(imds, diffIndexes, currIndexes);
            catch e
                throw(e)
            end
        end
        % Set method for ReadSize
        function set.ReadSize(imds, readSize)
            try
                iValidateReadSize(readSize);
                imds.ReadSize = readSize;
            catch e
                throw(e);
            end
        end
        % Set method for Labels
        function set.Labels(imds, labels)
            try
                labels = validateLabels(imds, labels);
                imds.Labels = labels;
            catch e
                throw(e);
            end
        end
        % Set CachedRead
        function set.CachedRead(imds, cachedRead)
            try
                cachedRead = validatestring(cachedRead, {'on', 'off'});
                imds.CachedRead = cachedRead;
            catch e
                throw(e);
            end
        end
        % Set MaxThreads
        function set.MaxThreads(imds, maxThreads)
            try
                maxThreads = iValidateMaxThreads(maxThreads);
                imds.MaxThreads = maxThreads;
                imds.IsMaxThreadsDefault = isequal(maxThreads, 'default');
            catch e
                rethrow(e);
            end
        end
        % Get Files
        function files = get.Files(imds)
            files = imds.Splitter.Files;
        end
        % Get Labels
        function labels = get.Labels(imds)
            labels = imds.Labels;
        end
        
        function set.PrefetchSize(imds, prefetchSize)            
            imds.PrefetchSize = prefetchSize;
        end
    end

    methods (Hidden)
         function imds = subsasgn(imds, S, B)
            % Store the indexes of added and emptied files.
            %   - AddedIndexes are used to add default labels.
            %   - EmptiedIndexes are used to empty out labels.
            try
                switch numel(S)
                    case 2
                        if isequal(S(1).type, '.')
                            switch S(1).subs
                                case 'Files'
                                    if isempty(B)
                                        imds.EmptyIndexes = S(2).subs{1};
                                    else
                                        idxes = S(2).subs{1};
                                        imds.AddedIndexes = idxes(idxes > imds.NumFiles);
                                    end
                            end
                        end
                    case 1
                        if isequal(S.type, '.')
                            switch S.subs
                                case 'Files'
                                    imds.FilesAssigned = true;
                            end
                        end
                end

                imds = builtin('subsasgn', imds, S, B);
            catch e
                throw(e)
            end
        end

        function initFromFileSplit(ds, filename, offset, len)
            files = ds.Files;
            l = ds.Labels;
            initFromFileSplit@matlab.io.datastore.CustomReadDatastore(ds, filename, offset, len);
            if ~isempty(l)
                newFiles = ds.Files;
                [~, index] = ismember(newFiles, files);
                if any(index == 0)
                    error(message('MATLAB:datastoreio:imagedatastore:invalidFilenameFromSplit', filename));
                else
                    ds.Labels = l(index);
                    setExistingCategories(ds);
                end
            end
            % set datastore not to use prefetch reading
            ds.IsReadFcnDefault = false;
        end
        
        function avgFs = getAverageFileSize(imds)
            avgFs = mean([imds.Splitter.Splits(:).FileSize]);            
        end
    end

    methods (Access = private)
        
        % Set IsReadFcnDefault if ReadFcn is the default readDatastoreImage private method
        function setIsReadFcnDefault(imds, readFcn)
            fcnInfo = functions(readFcn);
            pvtFile = fullfile(fileparts(mfilename('fullpath')), 'private', 'readDatastoreImage.m');
            tf = isfield(fcnInfo, 'file') && isequal(fcnInfo.file, pvtFile);
            imds.IsReadFcnDefault = tf && isfield(fcnInfo, 'parentage') && isequal(fcnInfo.parentage, {'readDatastoreImage'});
        end

        % Set File Extensions of any new files or remove file extensions
        % of files that are removed with subsasgn
        function setFileExtensions(imds, diffIndexes, currIndexes)
            fileExts = cell(imds.NumFiles, 1);
            if nargin == 1
                idxes = 1:imds.NumFiles;
            else
                fileExts(~diffIndexes) = imds.FileExts(currIndexes);
                idxes = find(diffIndexes);
            end
            for ii = 1:numel(idxes)
                filename = imds.Splitter.Splits(idxes(ii)).Filename;
                [~,~,ext] = fileparts(filename);
                fileExts{idxes(ii)} = lower(ext);
            end
            imds.FileExts = fileExts;
        end

        function [data, info] = preFetchRead(imds)
            if ~hasdata(imds)
                error(message('MATLAB:datastoreio:splittabledatastore:noMoreData'));
            end
            readSize = getTrueReadSize(imds);
            [filesToRead, idxes] = nextFilesToRead(imds, readSize);
            [data, imds.BatchReader] = readUsingPreFetcher(imds, filesToRead, idxes, imds.BatchReader, readSize);
            info = getInfoForBatch(imds, filesToRead, idxes);

            % Set the SplitIdx appropriately
            newSplitIdx = imds.SplitIdx + readSize;
            if newSplitIdx > imds.NumFiles
                imds.SplitIdx = imds.NumFiles;
                % For hasdata to return false
                imds.SplitReader.ReadingDone = true;
                return;
            end
            imds.SplitIdx = newSplitIdx;

            % start next set of files to read when the ImageDatastore
            % buffer is empty and prefetching is disabled
            if isempty(imds.DataBuffer) && imds.PrefetchSize <= imds.ReadSize
                [filesToRead, idxes] = nextFilesToRead(imds);
                formats = imds.FileExts(idxes);
                startBatchReading(imds, imds.BatchReader, filesToRead, formats);
            end
        end

        function [data, info] = readUsingSplitReader(imds, splitIndex)
            splitReader = copy(imds.SplitReader);
            splitReader.Split = imds.Splitter.Splits(splitIndex);
            reset(splitReader);
            [data, info] = getNext(splitReader);
            %Weizong modify here: part 1
            if (size(data,1)==227 || size(data,1)==224) && size(data,3)==1
                disp('WARNING-Weizong: auto convert image from gray to color')
                data=cat(3,data,data,data);
            end
            info.Label = getLabelUsingIndex(imds, splitIndex);
        end

        function [data, reader] = readUsingPreFetcher(imds, filesToRead, idxes, reader, readSize)
            if isempty(reader)
                reader = matlab.io.datastore.internal.BatchImreader;
                formats = imds.FileExts(idxes);
                if imds.PrefetchSize > imds.ReadSize
                    % Remote location and first prefetching. We need to
                    % bootstrap the prefetching process. The files
                    % downloaded from startBatchReading will be also read
                    % subsequently by the read function.
                    [filesToPrefetch, prefetchIdxes] = nextFilesToRead(imds, imds.PrefetchSize);
                    prefetchFormats = imds.FileExts(prefetchIdxes);
                    startBatchReading(imds, reader, filesToPrefetch, prefetchFormats);
                else
                    startBatchReading(imds, reader, filesToRead, formats);
                end
                
            end
            if isempty(imds.DataBuffer)
                % Read from the imds buffer
                [tmp_cell, imds.ErrorBuffer] = read(reader);
                %Weizong modify here: part 2
                if ~isempty(tmp_cell)
                    tmp_size=size(tmp_cell{1});
                    if tmp_size(1)==227 || tmp_size(1)==224
                        for iii=1:length(tmp_cell)
                            rand_rot=rand;
                            tmp_img=tmp_cell{iii};
                            if rand_rot<0.25
%                                 deg=0;
                            elseif rand_rot<0.5
%                                 deg=90;
                                tmp_img=rot90(tmp_img,1); %rot90 is much faster than imrotate
                            elseif rand_rot<0.75
%                                 deg=180;
                                tmp_img=rot90(tmp_img,2);
                            else
%                                 deg=270;
                                tmp_img=rot90(tmp_img,3);
                            end
                            tmp_img=tmp_img*(min(1,1-(rand-0.5)*0.5));%half chance it is 1, half chance within 0.75-1
                            if tmp_size(end)==tmp_size(1)
                                tmp_img=cat(3,tmp_img,tmp_img,tmp_img);
                            end
                            tmp_cell{iii}=tmp_img;
                        end
                        imds.DataBuffer=tmp_cell;
                    end
                end
                if imds.PrefetchSize > imds.ReadSize
                    % If we are prefetching, start a new prefetching
                    % phase.The next imds.read will be getting data from
                    % imds.DataBuffer, while startBatchReading is
                    % downloading the next images from the remote location
                    [filesToPrefetch, prefetchIdxes] = nextFilesToRead(imds, imds.PrefetchSize);
                    prefetchFormats = imds.FileExts(prefetchIdxes);
                    startBatchReading(imds, reader, filesToPrefetch, prefetchFormats);
                end
            end
            
            bufferSize = numel(imds.DataBuffer);
            if bufferSize < readSize
                remIdxes = bufferSize+1:readSize;
                formats = imds.FileExts(idxes(remIdxes));
                startBatchReading(imds, reader, filesToRead(remIdxes), formats);
                [d, e] = read(reader);
                %Weizong modify part 3 (not very important, but in case)
                if (size(d,1)==227 || size(d,1)==224) && size(d,3)==1
                    disp('WARNING-Weizong: bufferSize < readSize & convert gray to rgb')
                    d=cat(3,d,d,d);
                end
                data = [imds.DataBuffer; d];
                errors = [imds.ErrorBuffer; e];
                imds.DataBuffer = [];
                imds.ErrorBuffer = [];
            else
                data = imds.DataBuffer(1:readSize);
                errors = imds.ErrorBuffer(1:readSize);
                imds.DataBuffer(1:readSize) = [];
                imds.ErrorBuffer(1:readSize) = [];
            end
            errIdxes = find(errors);
            % When this is invoked from datastore/readall, idxes is a
            % logical vector. When invoked by datastore/read, this is
            % a double vector of the form start:end
            isReadAll = islogical(idxes);
            if ~isempty(errIdxes)
                for ii = 1:numel(errIdxes)
                    errIdx = errIdxes(ii);
                    try
                        % If reading all the files, ignore the split index,
                        % otherwise use it
                        if isReadAll
                            splitIdx = 0;
                        else
                            splitIdx = imds.SplitIdx - 1;
                        end
                        % readUsingSplitReader reads from remote locations
                        % as well. SplitIdx points to the current file to
                        % be read errIdxes are 1:n, n - number of files to
                        % be read
                        data{errIdx} = readUsingSplitReader(imds, splitIdx + errIdx);
                    catch ME
                        % Once errored, make BatchReader empty so we try reading next time.
                        imds.BatchReader = [];
                        msg = message('MATLAB:datastoreio:imagedatastore:unableToReadFile', filesToRead{errIdx});
                        mexc = MException(msg);
                        mexc = addCause(mexc, ME);
                        throw(mexc);
                    end
                end
            end
        end

        function startBatchReading(imds, reader, filesToRead, formats)
            if ~imds.IsMaxThreadsDefault
                reader.MaxThreads = imds.MaxThreads;
            end
            startRead(reader, filesToRead, formats);
        end

        function info = getInfoForBatch(imds, filesToRead, idxes)
            if nargin < 3
                [filesToRead, idxes] = nextFilesToRead(imds);
            end
            if imds.ReadSize == 1 && iscell(filesToRead)
                filesToRead = filesToRead{1};
            end
            info.Filename = filesToRead;
            fileSizes = [imds.Splitter.Splits(idxes).FileSize];
            info.FileSize = fileSizes(:);
            info.Label = getLabelUsingIndex(imds, idxes);
        end

        function readSize = getTrueReadSize(imds)
            readSize = imds.ReadSize;
            remSize = imds.NumFiles - imds.SplitIdx + 1;
            if readSize > remSize
                readSize = remSize;
            end
        end
        
        function idxes = nextIndexes(imds, readSize)
            % Return the next set of indices corresponding to the next batch of data
            
            if readSize > imds.ReadSize && ~isempty(imds.DataBuffer)
                % If prefetching, read file in advance, i.e., start =
                % current index + readSize
                sIdx = imds.SplitIdx + readSize;
            else
                % Otherwise return the current index
                sIdx = imds.SplitIdx;
            end
            eIdx = min(sIdx + readSize - 1, imds.NumFiles);
            idxes = sIdx:eIdx;
        end
    
        function [files, idxes] = nextFilesToRead(imds, readSize)
            % Return the next set of file and indices for the next batch of
            % data to read
            if nargin == 1
                readSize = getTrueReadSize(imds);
            end
            idxes = nextIndexes(imds, readSize);
            files = {imds.Splitter.Splits(idxes).Filename};
            files = files(:);
        end

        function initDatastore(imds, files, nv)
            import matlab.io.datastore.ImageDatastore;
            import matlab.io.datastore.internal.validators.validateCustomReadFcn;
            import matlab.io.datastore.internal.isIRI;

            validateCustomReadFcn(nv.ReadFcn, true, ImageDatastore.CONVENIENCE_CONSTRUCTOR_FCN_NAME);
            [~, files, fileSizes, imds.FileExts] = ImageDatastore.supportsLocation(files, nv);
            imds.FileExts = lower(imds.FileExts);
            imds.SplitterName = ImageDatastore.WHOLE_FILE_CUSTOM_READ_SPLITTER_NAME;

            imds.IsReadFcnDefault = ismember('ReadFcn', nv.UsingDefaults);
            
            if ~isempty(files) && any(isIRI(files) & ~isSupportedIRI(files))
                imds.IsReadFcnDefault = false;
            end
            
            imds.NumFiles = numel(files);
            if ~ismember('ReadSize', nv.UsingDefaults)
                iValidateReadSize(nv.ReadSize);
            end
            imds.ReadSize = nv.ReadSize;
            if ismember('Labels', nv.UsingDefaults) && strcmpi(nv.LabelSource, 'foldernames') && imds.NumFiles ~= 0
                % When 'LabelSource' is 'foldernames', Labels is a categorical.
                l = cellfun(@iLabelFromFile, files, 'UniformOutput', false);
                nv.Labels = categorical(l);
            end

            % initReadFcnAndFiles sets the datastore's ReadFcn and passes rest
            % of the varargin inputs to the splitter.
            nv.Files = files;
            nv.FileSizes = fileSizes;
            initReadFcnAndFiles(imds, [1:numel(files)], nv);
            imds.Labels = nv.Labels;
        end

        function resolveLabelsIndexes(imds)
            % If the filenames are modified, Labels are not changed
            % If any additional files are added, we add default values
            % for the respective Labels.
            %
            %    Labels Type   Default
            %    -----------   -------
            %    cellstr        ''
            %    numerical      0
            %    logical        false
            %    categorical    <undefined>
            %
            if isempty(imds.Labels)
                return;
            end
            if isempty(imds.Files)
                imds.Labels = {};
                return;
            end
            if ~isempty(imds.EmptyIndexes)
                imds.Labels(imds.EmptyIndexes) = [];
                setExistingCategories(imds);
            end
            if ~isempty(imds.AddedIndexes)
                nv = numel(imds.AddedIndexes);
                addLabels = iGetDefaultLabels(imds.Labels, nv);
                imds.Labels(imds.AddedIndexes) = addLabels;
            end
            if imds.FilesAssigned
                imds.Labels = {};
            end
        end

        function label = getLabelUsingIndex(imds, idx)
            if isempty(imds.Labels)
                label = {};
                return;
            end

            if iscell(imds.Labels)
                label = imds.Labels{idx};
            else
                label = imds.Labels(idx);
            end
        end

        function labels = validateLabels(imds, labels)
            if iscell(labels)
                if isempty(labels)
                    return;
                end
                if ~iscellstr(labels)
                    error(message('MATLAB:datastoreio:imagedatastore:labelsNotCellstr'));
                end
            end
            classes = {'numeric', 'cell', 'categorical', 'logical'};
            attrs = {'numel', imds.NumFiles};
            import matlab.io.datastore.ImageDatastore;
            validateattributes(labels, classes, attrs, ...
                ImageDatastore.CONVENIENCE_CONSTRUCTOR_FCN_NAME, 'Labels');
            labels = labels(:);
        end

        function setExistingCategories(imds)
            % Create categorical with only the existing categories, so that
            % the resulting labels will not have categories of the
            % Labels.
            l = imds.Labels;
            if isa(l, 'categorical')
                u = unique(l);
                imds.Labels = categorical(l, u(~isundefined(u)));
            end
        end

        function initializeSubsAsgnIndexes(imds)
            % subsasgn for Files is called before set.Files.
            % reset EmptyIndexes and AddedIndexes if there's an error,
            % so we don't change the labels.
            imds.EmptyIndexes = [];
            imds.AddedIndexes = [];
            imds.FilesAssigned = false;
        end

        function [u, c, g, n] = groupAndCountLabels(imds)
            l = imds.Labels;
            [g, u] = findgroups(l);
            c = splitapply(@numel, l, g);
            n = [];
            switch class(l)
                case 'categorical'
                    % count the number of undefined labels
                    [c, u, n] = iAddUngroupedCounts(isundefined(l), categorical(nan), c, u, n);
                case 'cell'
                    % count the number of empty strings
                    % empty string indexes are NaNs in grouping indexes
                    [c, u, n] = iAddUngroupedCounts(isnan(g), {''}, c, u, n);
                case 'logical'
                    % nothing to do for logical; all logicals are always grouped
                otherwise
                    % It has to be numerical at this point, since
                    % labels can only be numerical, logical, cellstr or categorical.
                    [c, u, n] = iAddUngroupedCounts(isnan(l), nan, c, u, n);
            end
        end

        function initReadFcnAndFiles(imds, indexes, info)
        %INITREADFCNANDFILES Initialize datastore with specific file indexes and ReadFcn.
        %   This can be used to initialize the datastore with ReadFcn and files/fileSizes
        %   found previously or already existing in the splitter information.
        %
        %   See also ReadFcn.
            narginchk(2,3);
            if nargin == 3
                readFcn = info.ReadFcn;
                files = info.Files;
                fileSizes = info.FileSizes;
                includeSubfolders = info.IncludeSubfolders;
            else
                % use splitter information from the info struct and initialize.
                readFcn = imds.Splitter.ReadFcn;
                if isempty(imds.Splitter.Splits)
                    fileSizes = [];
                    files = {};
                else
                    files = imds.Splitter.Files;
                    fileSizes = [imds.Splitter.Splits.FileSize];
                end
                % No need to lookup subfolders, if we already have files from the
                % splitter which are resolved.
                includeSubfolders = false;
            end

            if ~isempty(fileSizes)
                fileSizes = fileSizes(indexes);
            end

            if ~isempty(files)
                files = files(indexes);
            end
            % Signature for initFromReadFcn:
            %    initFromReadFcn(ds, readFcn, files, fileSizes, includeSubfolders)
            %
            %  - files, fileSizes - Pass files and filesizes to initialization of the splitter,
            %    so we don't lookup the path and verify the existence of files
            %  - includeSubfolders - true/false to include recursive sub-folders
            initFromReadFcn(imds, readFcn, files, fileSizes, includeSubfolders);
        end
    end

    methods (Access = protected)

        function validateReadFcn(imds, readFcn)

            % validateReadFcn is called from set.ReadFcn
            import matlab.io.datastore.ImageDatastore;
            import matlab.io.datastore.internal.validators.validateCustomReadFcn;
            validateCustomReadFcn(readFcn, false, ImageDatastore.CONVENIENCE_CONSTRUCTOR_FCN_NAME);

            % Set the private IsReadFcnDefault value
            setIsReadFcnDefault(imds, readFcn);
        end

        function displayScalarObject(imds)
            % header
            disp(getHeader(imds));
            group = getPropertyGroups(imds);
            filesEmpty = imds.NumFiles == 0;
            labels = imds.Labels;
            labelsEmpty = isempty(labels);
            if ~filesEmpty || ~labelsEmpty
                detailsStr = evalc('details(imds)');
                nsplits = strsplit(detailsStr, '\n');
            end
            import matlab.io.datastore.internal.cellArrayDisp;
            import matlab.io.datastore.internal.vectorDisp;
            if ~filesEmpty
                filesStr = nsplits(~cellfun(@isempty, strfind(nsplits, 'Files: ')));
                % Find the indent spaces from details
                nFilesIndent = strfind(filesStr{1}, 'Files: ') - 1;
                if nFilesIndent > 0
                    % File Properties
                    filesIndent = [sprintf(repmat(' ',1,nFilesIndent)) 'Files: '];
                    nlspacing = sprintf(repmat(' ',1,numel(filesIndent)));
                    filesStrDisp = cellArrayDisp(imds.Files, true, nlspacing);
                    disp([filesIndent filesStrDisp]);
                    % Remove Files property from the group, since custom
                    % display is used for Files.
                    group.PropertyList = rmfield(group.PropertyList, 'Files');
                end
            end
            if ~labelsEmpty && imds.NumFiles > 1
                labelsStr = nsplits(~cellfun(@isempty, strfind(nsplits, 'Labels: ')));
                nLabelsIndent = strfind(labelsStr{1}, 'Labels: ') - 1;
                if nLabelsIndent > 0
                    labelsIndent = [sprintf(repmat(' ',1,nLabelsIndent)) 'Labels: '];
                end
                if iscell(labels)
                    labelsStrDisp = cellArrayDisp(labels, false, '');
                else
                    labelsStrDisp = vectorDisp(labels);
                end
                disp([labelsIndent, labelsStrDisp]);
                group.PropertyList = rmfield(group.PropertyList, 'Labels');
            end
            matlab.mixin.CustomDisplay.displayPropertyGroups(imds, group);
            disp(getFooter(imds));
        end
    end

    methods (Hidden)
        function files = getUnresolvedFiles(imds)
            files = imds.UnResolvedFiles;
        end
    end

    methods (Static, Hidden)

        function varargout = supportsLocation(loc, nvStruct)
            % This function is responsible for determining whether a given
            % location is supported by ImageDatastore. It also returns a
            % resolved filelist.
            i = imformats;
            defaultExtensions = strcat('.', [i.ext]);
            [varargout{1:nargout}] = matlab.io.datastore.FileBasedDatastore.supportsLocation(loc, nvStruct, defaultExtensions);
        end

        function outds = loadobj(ds)
            % ImageDatastore was introduced in 15b
            originatingVersion = '2015b';
            currVersion = version('-release');
            switch class(ds)
                case 'struct'
                    % load datastore from struct
                    if isfield(ds, 'SchemaVersion') && ~isempty(ds.SchemaVersion)
                        originatingVersion = ds.SchemaVersion;
                    end
                    ds = matlab.io.datastore.ImageDatastore.loadFromStruct(ds);
                case 'matlab.io.datastore.ImageDatastore'
                    if isprop(ds, 'SchemaVersion') && ~isempty(ds.SchemaVersion)
                        originatingVersion = ds.SchemaVersion;
                    end
            end
            switch originatingVersion
                case {'2015b', '2016a'}
                    % ReadSize was introduced in 16b
                    ds.ReadSize = 1;
                    if isequal(originatingVersion, '2015b') && isequal(currVersion, '2016a')
                        % 15b in 16a: Labels was introduced in 16a
                        ds.Labels = {};
                    end
            end
            fcnInfo = functions(ds.ReadFcn);
            % if saved version is less than 2017a,
            if ~iCompareVersion(originatingVersion, '2017a')
                % Re-assign default ReadFcn from ImageDatastore's
                % private readDatastoreImage
                if isfield(fcnInfo, 'class') && isequal(fcnInfo.class,  'matlab.io.datastore.ImageDatastore') ...
                    && isequal(fcnInfo.function, 'readDatastoreImage')
                    ds.ReadFcn = @readDatastoreImage;
                end
            end

            if isempty(ds.NumFiles)
                ds.NumFiles = numel(ds.Files);
            end
            if isempty(ds.FileExts)
                setFileExtensions(ds);
            end

            % Set if ReadFcn is the default in the saved datastore
            setIsReadFcnDefault(ds, ds.ReadFcn);

            if ds.Splitter.NumSplits ~= 0
                % create a split reader that points to the
                % first split index.
                if ds.SplitIdx == 0
                    ds.SplitIdx = 1;
                end
                % create a stub reader so copy() works fine as it expects
                % a non empty datastore to have a reader.
                ds.SplitReader = ds.Splitter.createReader(ds.SplitIdx);
            end
            outds = loadobj@matlab.io.datastore.FileBasedDatastore(ds);
        end
    end

    methods (Static, Access = private)
        function ds = loadFromStruct(dsStruct)

            % empty datastore
            ds = matlab.io.datastore.ImageDatastore({});
            ds.Splitter = dsStruct.Splitter;
            initReadFcnAndFiles(ds, [1:ds.Splitter.NumSplits]);

            ds.BatchReader = [];
            ds.DataBuffer = [];
            ds.ErrorBuffer = [];
            fieldsToRemove = {'Splitter', 'BatchReader',...
                              'DataBuffer', 'ErrorBuffer',...
                              };
            fieldList = fields(dsStruct);
            % Previous to 16b, BatchReader, DataBuffer, ErrorBuffer might not be present
            fieldsToRemove = iIntersectStrings(fieldList, fieldsToRemove);

            if ~isempty(fieldsToRemove)
                dsStruct = rmfield(dsStruct, fieldsToRemove);
                fieldList = fields(dsStruct);
            end

            for fieldIndex = 1: length(fieldList)
                field = fieldList{fieldIndex};
                ds.(field) = dsStruct.(field);
            end
        end
    end
end

function parsedStruct = iParseNameValues(varargin)
    persistent inpP;
    if isempty(inpP)
        import matlab.io.datastore.ImageDatastore;
        inpP = inputParser;
        addParameter(inpP, 'ReadSize', 1);
        addParameter(inpP, 'ReadFcn', @readDatastoreImage);
        addParameter(inpP, 'Labels', {});
        addParameter(inpP, 'LabelSource', 'none', @(x)validateattributes(x, {'char'}, {'nonempty'}));
        addParameter(inpP, 'IncludeSubfolders', false);
        addParameter(inpP, 'FileExtensions', -1);
        inpP.FunctionName = ImageDatastore.CONVENIENCE_CONSTRUCTOR_FCN_NAME;
    end
    parse(inpP, varargin{:});
    parsedStruct = inpP.Results;
    parsedStruct.LabelSource = validatestring(parsedStruct.LabelSource, {'none', 'foldernames'});
    parsedStruct.UsingDefaults = inpP.UsingDefaults;
    if ~ismember('Labels', inpP.UsingDefaults) && ~ismember('LabelSource', inpP.UsingDefaults) ...
        && strcmpi(parsedStruct.LabelSource, 'foldernames')
        error(message('MATLAB:datastoreio:imagedatastore:labelsLabelSourceCombined'));
    end
end

function label = iLabelFromFile(file)
    [~, label] = fileparts(fileparts(file));
end

function labels = iGetDefaultLabels(sampleLabels, numLabels)
    switch class(sampleLabels)
        case 'cell'
            labels = repmat({''}, numLabels, 1);
        case 'logical'
            labels = false(numLabels,1);
        case 'categorical'
            labels = categorical(nan(numLabels,1));
        otherwise
            % It has to be numerical at this point, since
            % labels can only be numerical, logical, cellstr or
            % categorical.
            labels = zeros(numLabels, 1,'like',sampleLabels);
    end
end

% Add ungrouped counts and values to the count c, unique groups u, and
% the ungrouped indexes to n;
function [c, u, n] = iAddUngroupedCounts(ungrouped, ungroupedValue, c, u, n)
    ungroupedCount = nnz(ungrouped);
    if ungroupedCount ~= 0
        c = [c; ungroupedCount];
        u = [u; ungroupedValue];
        n = ungrouped;
    end
end

function maxThreads = iValidateMaxThreads(maxThreads)
    if ischar(maxThreads)
        maxThreads = validatestring(maxThreads, {'default'}, mfilename, 'MaxThreads');
        return;
    end
    classes = {'numeric'};
    attrs = {'scalar', 'positive', 'integer'};
    import matlab.io.datastore.ImageDatastore;
    validateattributes(maxThreads, classes, attrs, mfilename, 'MaxThreads');
end

function iValidateReadSize(rsize)
    classes = {'numeric'};
    attrs = {'>=', 1, 'scalar', 'positive', 'integer', 'nonsparse'};
    import matlab.io.datastore.ImageDatastore;
    validateattributes(rsize, classes, attrs, ...
        ImageDatastore.CONVENIENCE_CONSTRUCTOR_FCN_NAME, 'ReadSize');
end

% Compares MATLAB versions vOne and vTwo
% vOne, vTwo are versions obtained from version('-release') command.
% Returns true iff vOne >= vTwo
function tf = iCompareVersion(vOne, vTwo)
tf = true;
if isequal(vOne, vTwo)
    return;
end
vOneNum = str2num(vOne(1:4));
vTwoNum = str2num(vTwo(1:4));
if vOneNum > vTwoNum
    return;
elseif vOneNum == vTwoNum
    tf = vOne(5) > vTwo(5);
    return;
end
tf = false;
end

% A for loop version of intersect. Remove string items from setTwo
% if not present in setOne argument.
function setTwo = iIntersectStrings(setOne, setTwo)
    num = numel(setTwo);
    idxes = false(num, 1);
    for ii = 1:num
        idxes(ii) = any(strcmp(setOne, setTwo(ii)));
    end
    setTwo(~idxes) = [];
end
