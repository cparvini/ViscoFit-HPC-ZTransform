function [dataStruct] = LoadAFMData(pathname,varargin)
%LOADAFMDATA Load AFM Data from an External Source
%   This function takes in a pathname, and will then parse it looking for
%   relevant AFM SFS files. Implementing new filetypes involves adding new
%   cases to the switch-case that ties the label (i.e. first string in the
%   filename, which is underscore "_" separated) to the correct steps for
%   processing. For example, if you are attempting to load
%   AFMData_Experiment1.ibw, then there must be a case for when strSwitch
%   is "AFMData". That case should (for the first switch-case) determine
%   which files in the directory are the correct ones to load. Then, the
%   second switch case requires the user to write the code required to load
%   at least the following variables: z (Z-Sensor [m]), d (Deflection [m]),
%   r_tip (Tip Size, [m]; Note: this can also be cone angle [deg] for cone
%   tips, but for simplicity the variable should still be named "r_tip"),
%   nu_sample (Poisson's ratio for the sample if known, if not set to 0.5
%   for incompressible), k_cantilever (stiffness of the cantilever [N/m]),
%   t (time array [s]), dt (timestep [s]), and N_sam (Number of samples in
%   the t, z, and d arrays, used for bookkeeping). These variables must be
%   stored in rows (one for each file) in a structure named "dataStruct".
%   The file will then organize and clean up the input data, find the
%   region of interest for fitting, and also perform averaging operations
%   for all of the files that roughly share an approach velocity. The
%   original data is kept in the first 1-N rows (where N is the number of
%   files provided), and the new "averaged" datasets are stored at the END
%   of dataStruct. There is one row for each approach velocity detected. 

% Default values
includeRetract = 0;             % Include data from the retract curve
filterType = 'none';            % Choose the filter used to smooth data
N = 2;                          % Order of Butterworth filter, if used
cutoff_Hz = 5000;               % Cutoff frequency
findRep = 'forward';            % Search direction for the repulsive region
removeNegatives = true;         % Remove negative values in the data stream
createAverage = true;           % Create the averaged rows for each approach vel.

% Read varargin values
if ~isempty(varargin)
    if ~isa(varargin{1},'struct')
        error('You can only pass a STRUCT containing settings to LoadAFMData(). Please verify you have done so, because the type does not appear to be "struct".');
    else
        % Load the settings
        inputSettings = varargin{1};
        
        % Load mandatory settings
        try
            includeRetract = inputSettings.includeRetract;
        end
        try
            filterType = inputSettings.filterType;
        end
        try
            findRep = inputSettings.findRep;
        end
        try
            removeNegatives = inputSettings.removeNegatives;
        end
        try
            createAverage = inputSettings.createAverage;
        end
        
        % Load conditional settings
        if strcmpi(filterType,'butter')
            N = inputSettings.N;
            cutoff_Hz = inputSettings.cutoff_Hz;
        end
    end
end

FilesCheck=dir([pathname filesep '*.*']);

% Remove Directories
FilesCheck=FilesCheck(~ismember({FilesCheck.name},{'.','..'}));
toRemove = find([FilesCheck.isdir] == 1);
FilesCheck(toRemove) = [];

% Remove Filetypes To Ignore
toRemove = find(~endsWith({FilesCheck.name}, {'.ibw','.txt','.spm','.mat','.csv','.jpk-qi-data'}));
FilesCheck(toRemove) = [];

toRemove = find(contains({FilesCheck.name}, {'settingsStruct','Settings'}));
FilesCheck(toRemove) = [];

toRemove = find(contains({FilesCheck.name}, {'FitResults','PlotResults','mapStruct','log.txt'}));
FilesCheck(toRemove) = [];

for i = 1:length(FilesCheck)
    FilesTempPrep = FilesCheck(i).name;
    FilesTempPrep = strsplit(FilesTempPrep, {'_' '.' ' '},'CollapseDelimiters',true);
    strSwitchPrep = FilesTempPrep{1};
    tempLabels{i} = strSwitchPrep;
end

if ~exist('tempLabels','var')
    error('ERROR: No valid AFM files found! Please ensure you have selected the correct directory, and that LoadAFMData contains a case describing how to load that filetype.');
end

if length(FilesCheck) > 1
    FilesTemp = FilesCheck(1).name;
    FilesTemp = strsplit(FilesTemp, {'_' '.' ' '},'CollapseDelimiters',true);
    
    strSwitch = FilesTemp{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            Files=dir([pathname '/*.ibw']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
        
        case lower('TestCondition')
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
            Files(toRemove) = [];
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(FileInfo{5})*1E-9;                  % Approach Velocity, nm/s
                point_number(k) = k;                                        % Position Number
                run_number(k) = 1;                                          % Run Number
            end
            
        case lower({'HPAF','HPDE'})
            Files = dir([pathname '/*.jpk-qi-data']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
            Files(toRemove) = [];
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);
            end
            
        otherwise
            if endsWith(FilesCheck.name, {'.jpk-qi-data'})
                Files = dir([pathname '/*.jpk-qi-data']);
                toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
                Files(toRemove) = [];
                for k=1:length(Files)
                    FileNames = Files(k).name;
                    FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);
                end
            else
               error('There was no built-in case for the files you have in the chosen directory. Make sure you have added a case to LoadAFMData() for your unique identifier, or are analyzing QI maps.'); 
            end
        
    end
    
else
    
    FilesTemp = FilesCheck.name;
    FilesTemp = strsplit(FilesTemp, {'_' '.' ' '},'CollapseDelimiters',true);
    
    strSwitch = FilesTemp{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            Files=dir([pathname '/*.ibw']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number

        case lower('TestCondition')
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
            Files(toRemove) = [];
            
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

            v_approach = str2num(FileInfo{5})*1E-9;                 % Approach Velocity, nm/s
            point_number = 1;                                       % Position Number
            run_number = 1;                                         % Run Number
            
        case lower({'HPAF','HPDE'})
            Files = dir([pathname '/*.jpk-qi-data']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
            Files(toRemove) = [];
            
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

        otherwise
            if endsWith(FilesCheck.name, {'.jpk-qi-data'})
                Files = dir([pathname '/*.jpk-qi-data']);
                toRemove = find(contains({Files.name}, {'settingsStruct','Settings','FitResults'}));
                Files(toRemove) = [];

                FileNames = Files.name;
                FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);
            else
               error('There was no built-in case for handling the files you have in the chosen directory. Make sure you have added a case to LoadAFMData() for your unique identifier, or are analyzing QI maps.'); 
            end
            
    end
end

% Initialize our output structure
dataStruct = struct;

% Remove the files we don't care about
FilesRemove=(~ismember({Files.name},{FilesCheck.name}));
Files(FilesRemove) = [];

% Loop through the files and proceed with loading, parsing, and calculation
% for each one.
for k = 1:length(Files)
    FileNames = Files(k).name;
    FileInfo = strsplit(FileNames, {'_' '.' ' '},'CollapseDelimiters',true);
    
    strSwitch = FileInfo{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            % This case shows the steps used to load .ibw files from an
            % MFP3D running the Asylum AFM software. The file header is
            % read into headerValue and all of the relevant settings are
            % extracted from there. The actual data itself relies on the
            % function IBWread to load the AFM experiment observables. Note
            % that IBWread() was NOT created by the authors of this
            % repository/manuscript. A license for that script is included
            % in this directory, and named IBWread_license.txt. Credit goes
            % to Jakub Bialek (2009).
            RawData = IBWread([Files(k).folder filesep Files(k).name]);
            [headerValue,~] = strsplit(RawData.WaveNotes,'\r',...
            'DelimiterType','RegularExpression');

            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, m/s
            point_number_temp = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number_temp = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
            dataStruct(k).z = RawData.y(:,1);
            dataStruct(k).d = RawData.y(:,2);
            
            dataStruct(k).r_tip = input(sprintf('Please enter the tip radius used for the file "%s": ',Files(k).name));
            dataStruct(k).nu_sample = input(sprintf('Please enter the sample Poissons Ratio (nu) for file "%s": ',Files(k).name));
            
            % Find Spring Constant
            varIndex = find(contains(headerValue,'SpringConstant'),1);
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            k_cantilever(k) = str2num(temp{2}); % Value in N/m
            dataStruct(k).k_cantilever = k_cantilever(k);
            
            % Find Deflection InvOLS
            varIndex = find(contains(headerValue,'InvOLS'));
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            defl_InVOLS(k) = str2num(temp{2}); % Value in nm/V
            
            % Create time array
            varIndex = find(contains(headerValue,'NumPtsPerSec'));
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            dataStruct(k).dt = 1/str2num(temp{2}); % Value in s
            dataStruct(k).n_sam = RawData.Nsam;
            dataStruct(k).t = dataStruct(k).dt.*((1:size(dataStruct(k).d))-1)';
    
        case lower('TestCondition')
            % This case is for loading the "TestCondition" simulation files
            % which were generated for testing the code bases. The files
            % have a regular, relatively simple format and have two
            % separate .mat files for each: one contains the data from the
            % simulation (loaded here into RawData), and the other is a
            % settings file (used to perform the simulation, stored here in
            % settingsData).
            RawData = load([Files(k).folder filesep Files(k).name],'z','d','t','F');
                        
            settingsCheck=dir([pathname '/*.*']);

            % Remove Directories
            settingsCheck=settingsCheck(~ismember({settingsCheck.name},{'.','..'}));
            toRemove = find([settingsCheck.isdir] == 1);
            settingsCheck(toRemove) = [];

            % Remove Filetypes To Ignore
            toRemove = find(~endsWith({settingsCheck.name}, {'.mat'}));
            settingsCheck(toRemove) = [];

            toRemove = find(~contains({settingsCheck.name}, {'Settings'}));
            settingsCheck(toRemove) = [];
            
            if size(settingsCheck,1) > 1
                setInd = k;
            else
                setInd = 1;
            end
            
            settingsData = load([settingsCheck(setInd).folder filesep settingsCheck(setInd).name]);
            
            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;    % Approach Velocity, m/s
            point_number_temp = k;                                          % Position Number
            run_number_temp = 1;                                            % Run Number

            dataStruct(k).z = -(RawData.z-RawData.z(1));
            dataStruct(k).t = RawData.t;
            dataStruct(k).n_sam = numel(dataStruct(k).t);
            
            k_cantilever(k) = settingsData.settingsStruct.k_m1;
            dataStruct(k).k_cantilever = k_cantilever(k);
            dataStruct(k).r_tip = settingsData.settingsStruct.r_tip;
            dataStruct(k).nu_sample = settingsData.settingsStruct.nu_sample;
            dataStruct(k).d = RawData.F./k_cantilever(k);
            defl_InVOLS(k) = 1;
            dataStruct(k).dt = settingsData.settingsStruct.dt;
            
        case lower({'HPAF','HPDE'})
            % This is a test case for the map analysis functionality. The
            % data used here was collected using a JPK NanoWizard BioAFM
            % and a conical, live cell probe. This type of data loading is
            % unique from the legacy methods above, because there are a
            % large number of pixels being imported as opposed to
            % individual files corresponding to force curve experiments.
            % This is because the JPK Quantititative Imaging (QI) method
            % was used to collect both map and force spectroscopy data
            % simultaneously. This vastly increases the amount of data. 
            % Using the "open_JPK" function which has been adapted from the
            % original version on the MEX from Dr. Ortuso, R.D., each pixel
            % is introduced as a row in the dataStruct, with an ID column
            % that will store separate maps with unique IDs so the data can
            % be separated later on.
            oldMaps = dir([Files(k).folder sprintf('/mapStruct-%s-%s.mat',FileInfo{1},FileInfo{2})]);
            
            if isempty(oldMaps)
                tic
                fileStruct = open_JPK([Files(k).folder filesep Files(k).name]);
                loadTime = toc;
                fprintf('\nIt took %6.2f minutes to load %s.\n',loadTime/60,Files(k).name)
                save([Files(k).folder filesep sprintf('mapStruct-%s-%s.mat',FileInfo{1},FileInfo{2})],'fileStruct','-v7.3');
            else
                load([oldMaps.folder filesep oldMaps.name],'fileStruct'); % Load the previous filestruct!
            end
            numPixels = length(fileStruct);
            indShift = size(dataStruct,2);
            
            if indShift == 1
                indShift = 0;
            end
            
            for i_pix = (1:numPixels)
                % Store MapID
                dataStruct(i_pix+indShift).mapID = k;
                
                % Store the Z-Extension measurements from the map
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'z'),1);
                dataStruct(i_pix+indShift).z = fileStruct{i_pix}(idx).extend;
                dataStruct(i_pix+indShift).z = vertcat(dataStruct(i_pix+indShift).z,fileStruct{i_pix}(idx).retract);

                % Store the deflection
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'d'),1);
                dataStruct(i_pix+indShift).d = fileStruct{i_pix}(idx).extend;
                dataStruct(i_pix+indShift).d = vertcat(dataStruct(i_pix+indShift).d,fileStruct{i_pix}(idx).retract);
                
                % Store the Force
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'F'),1);
                dataStruct(i_pix+indShift).F = fileStruct{i_pix}(idx).extend;
                dataStruct(i_pix+indShift).F = vertcat(dataStruct(i_pix+indShift).F,fileStruct{i_pix}(idx).retract);

                % Store the Pixel Height
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'height'),1);
                dataStruct(i_pix+indShift).height = fileStruct{i_pix}(idx).extend;
                
                % Store the time array. Use the timestep (dt) to calculate
                % the retract time, too.
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'t'),1);
                dataStruct(i_pix+indShift).t = fileStruct{i_pix}(idx).extend;
                if isrow(dataStruct(i_pix+indShift).t)  dataStruct(i_pix+indShift).t =  dataStruct(i_pix+indShift).t'; end
                dataStruct(i_pix+indShift).dt = mode(round(gradient(dataStruct(i_pix+indShift).t),3,'significant'));
                dataStruct(i_pix+indShift).t = vertcat( dataStruct(i_pix+indShift).t,...
                    (dataStruct(i_pix+indShift).dt.*(1:length(fileStruct{i_pix}(idx).retract))'+dataStruct(i_pix+indShift).t(end)) );
                
                v_approach(i_pix+indShift) = round(mean(abs(gradient(dataStruct(i_pix+indShift).z)./dataStruct(i_pix+indShift).dt)),2,'significant');  % Approach Velocity, m/s
                
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'pixel'),1);
                point_number(i_pix+indShift) = str2double(fileStruct{i_pix}(idx).extend);     % Position Number
                run_number(i_pix+indShift) = 1;     % Run Number

                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'mapsize'),1);
                dataStruct(i_pix+indShift).mapSize = fileStruct{i_pix}(idx).extend;
                
                idx = find(contains({fileStruct{i_pix}(:).Channel_name},'stiffness'),1);
                k_cantilever(i_pix+indShift) = fileStruct{i_pix}(idx).extend;
                
                dataStruct(i_pix+indShift).r_tip = 15;   % Half-Angle, Degrees, of the Conical Tip
                dataStruct(i_pix+indShift).nu_sample = 0.5;  % Poisson's Ratio of the sample
                
            end
            
            clearvars fileStruct
            
        otherwise
            FileInfo = strsplit(FileNames, {'.'},'CollapseDelimiters',true);
            if endsWith(Files(k).name, {'.jpk-qi-data'})
                % This is the case for JPK QI mapping analysis. The data
                % used here was collected using a JPK NanoWizard BioAFM
                % and a conical, live cell probe. This type of data loading is
                % unique from the legacy methods above, because there are a
                % large number of pixels being imported as opposed to
                % individual files corresponding to force curve experiments.
                % This is because the JPK Quantititative Imaging (QI) method
                % was used to collect both height and force spectroscopy data
                % simultaneously. This vastly increases the amount of data. 
                % Using the "open_JPK" function, which has been adapted from the
                % original version on the MEX from Dr. Ortuso, R.D., each pixel
                % is introduced as a row in the dataStruct, with an ID column
                % that will store separate maps with unique IDs so the data can
                % be separated later on.
                oldMaps = dir([Files(k).folder sprintf('/mapStruct-%s.mat',FileInfo{1})]);

                if isempty(oldMaps)
                    tic
                    fileStruct = open_JPK([Files(k).folder filesep Files(k).name]);
                    loadTime = toc;
                    fprintf('\nIt took %6.2f minutes to load %s.\n',loadTime/60,Files(k).name)
                    save([Files(k).folder filesep sprintf('mapStruct-%s.mat',FileInfo{1})],'fileStruct');
                else
                    load([oldMaps.folder filesep oldMaps.name],'fileStruct'); % Load the previous filestruct!
                end
                numPixels = length(fileStruct);
                indShift = size(dataStruct,2);

                if indShift == 1
                    indShift = 0;
                end

                for i_pix = (1:numPixels)
                    % Store MapID
                    dataStruct(i_pix+indShift).mapID = k;

                    % Store the Z-Extension measurements from the map
                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'z'),1);
                    dataStruct(i_pix+indShift).z = fileStruct{i_pix}(idx).extend;
                    dataStruct(i_pix+indShift).z = vertcat(dataStruct(i_pix+indShift).z,fileStruct{i_pix}(idx).retract);

                    % Store the deflection
                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'d'),1);
                    dataStruct(i_pix+indShift).d = fileStruct{i_pix}(idx).extend;
                    dataStruct(i_pix+indShift).d = vertcat(dataStruct(i_pix+indShift).d,fileStruct{i_pix}(idx).retract);

                    % Store the Force
                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'F'),1);
                    dataStruct(i_pix+indShift).F = fileStruct{i_pix}(idx).extend;
                    dataStruct(i_pix+indShift).F = vertcat(dataStruct(i_pix+indShift).F,fileStruct{i_pix}(idx).retract);

                    % Store the Pixel Height
                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'height'),1);
                    dataStruct(i_pix+indShift).height = fileStruct{i_pix}(idx).extend;

                    % Store the time array. Use the timestep (dt) to calculate
                    % the retract time, too.
                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'t'),1);
                    dataStruct(i_pix+indShift).t = fileStruct{i_pix}(idx).extend;
                    if isrow(dataStruct(i_pix+indShift).t)  dataStruct(i_pix+indShift).t =  dataStruct(i_pix+indShift).t'; end
                    dataStruct(i_pix+indShift).dt = mode(round(gradient(dataStruct(i_pix+indShift).t),3,'significant'));
                    dataStruct(i_pix+indShift).t = vertcat( dataStruct(i_pix+indShift).t,...
                        (dataStruct(i_pix+indShift).dt.*(1:length(fileStruct{i_pix}(idx).retract))'+dataStruct(i_pix+indShift).t(end)) );

                    v_approach(i_pix+indShift) = round(mean(abs(gradient(dataStruct(i_pix+indShift).z)./dataStruct(i_pix+indShift).dt)),2,'significant');  % Approach Velocity, m/s

                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'pixel'),1);
                    point_number(i_pix+indShift) = str2double(fileStruct{i_pix}(idx).extend);     % Position Number
                    run_number(i_pix+indShift) = 1;     % Run Number

                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'mapsize'),1);
                    dataStruct(i_pix+indShift).mapSize = fileStruct{i_pix}(idx).extend;

                    idx = find(contains({fileStruct{i_pix}(:).Channel_name},'stiffness'),1);
                    k_cantilever(i_pix+indShift) = fileStruct{i_pix}(idx).extend;

                    dataStruct(i_pix+indShift).r_tip = 15;   % Half-Angle, Degrees, of the Conical Tip
                    dataStruct(i_pix+indShift).nu_sample = 0.5;  % Poisson's Ratio of the sample

                end

                clearvars fileStruct
                
            else
               error('There was no built-in case for handling the files you have in the chosen directory. Make sure you have added a case to LoadAFMData() for your unique identifier, or are analyzing QI maps.'); 
            end
                        
    end
    
end

for k = 1:size(dataStruct,2)

    if isfield(dataStruct,'mapID')
        numPixels = dataStruct(k).mapSize(1)*dataStruct(k).mapSize(2);
        ai = 1+numPixels*(dataStruct(k).mapID-1);
        bi = numPixels*(dataStruct(k).mapID);
    else
        ai = 1;
        bi = size(dataStruct,2);
    end
    
    % Find the approach portion of the data    
    [~, z_max_ind] = max(dataStruct(k).z);
    
    % Make sure we are handling row vectors
    if ~isrow(dataStruct(k).z) dataStruct(k).z = dataStruct(k).z'; end
    if ~isrow(dataStruct(k).d) dataStruct(k).d = dataStruct(k).d'; end
    if ~isrow(dataStruct(k).t) dataStruct(k).t = dataStruct(k).t'; end
    
    % Trim the data. If we include the retract portion of the dataset, then
    % we will go until the applied force detected is negative. If it is not
    % included, we end at the point of maximum z-sensor (maximum
    % deformation condition).
    if ~includeRetract
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    else
        F_temp = dataStruct(k).d(z_max_ind:end) .* k_cantilever(k);
        if ~isempty(F_temp) && length(F_temp) > 1
            non_contact_ind = find(F_temp < 0,1);
            if isempty(non_contact_ind) non_contact_ind = length(F_temp)-1; end
        else
            non_contact_ind = 0;
        end
        z_max_ind = non_contact_ind+z_max_ind;
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    end
    
    % Filter and Shift z_sensor data
    switch filterType
        case 'butter'
            % Create the butterworth
            [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(k).dt)),'low'); % This makes a lowpass filter

            d_smooth = (filter(b,a,dataStruct(k).d_approach));             % Next, apply the filter
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end
            
            % Calculate Deflection Offset
            [~, d_min_ind] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
        case {'FIR','IIR'}
            
            Fs = 1/(dataStruct(k).dt); 
            Fstop = ( (round((v_approach(k)),2,'significant')...
                /round(max(v_approach(ai:bi)),2,'significant'))...
                /dataStruct(k).dt )/10;
            if Fstop >= 1/(2*dataStruct(k).dt)
                Fstop = 1/(2.05*dataStruct(k).dt);
            elseif Fstop < 1/(10*dataStruct(k).dt)
                Fstop = 1/(10*dataStruct(k).dt);
            end
            Fpass = Fstop*0.01;
            Rp = 0.01;
            Astop = 80;
%             LPF = dsp.LowpassFilter('SampleRate',Fs, ...
%                                      'FilterType',filterType, ...
%                                      'PassbandFrequency',Fpass, ...
%                                      'StopbandFrequency',Fstop, ...
%                                      'PassbandRipple',Rp, ...
%                                      'StopbandAttenuation',Astop);
%             delay = floor(mean(grpdelay(LPF)));
%             d_smooth = LPF(dataStruct(k).d_approach);

            Nfilter = round((Fs/(Fstop-Fpass))*(Astop/22)); % "Fred Harris Rule of Thumb" for Filter Order
            firlo = fir1(Nfilter,0.5,'low');
            d_smooth = filter(firlo,1,dataStruct(k).d_approach);
            delay = floor(mean(grpdelay(firlo,Nfilter)));
            
            sf = d_smooth;
            if delay < length(sf)
                % Correct filter delay
                sf(1:delay) = [];
            else
                delay = 0;
            end
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0 || dSmoothMin > length(dataStruct(k).z_approach)
                dSmoothMin = 1;
            end

            d_smooth = sf;
            
            % Calculate Deflection Offset
            [~, d_min_ind] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected((delay+1):end);

        case 'none'
            d_smooth = dataStruct(k).d_approach;
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end
            
            % Calculate Deflection Offset
            [~, d_min_ind] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
    end
    
    % Correct the z-sensor data
    dataStruct(k).z_full_corrected = dataStruct(k).z - ...
        dataStruct(k).z_approach(dSmoothMin) + ...
        dataStruct(k).d_corrected(dSmoothMin);

    % Store our detected offsets
    dataStruct(k).d0 = dataStruct(k).d_corrected(dSmoothMin);
    dataStruct(k).z0 = dataStruct(k).z_corrected(dSmoothMin);
    
    % Calculate Force and Indentation
    dataStruct(k).F = dataStruct(k).d_corrected .* k_cantilever(k);
    dataStruct(k).F_smooth = d_smooth .* k_cantilever(k);
    dataStruct(k).d_smooth = d_smooth;
    dataStruct(k).z_smooth = z_smooth;
    dataStruct(k).h = (dataStruct(k).z(1:z_max_ind) - dataStruct(k).z0)...
        - (dataStruct(k).d(1:z_max_ind) - dataStruct(k).d0);
    
    % Get Repulsive Force Application Portion of the Data
    n_offset = length(dataStruct(k).d_corrected(dSmoothMin:z_max_ind));
    n_offset_smooth = length(dataStruct(k).d_corrected(dSmoothMin:(z_max_ind-delay)));

    if n_offset < 10
        fprintf('\nBad Pixel: %d\nStoring fake data here.\n',k);
        
        % Store fake arrays so the script doesn't fail for this bad pixel.
        dataStruct(k).t_r = dataStruct(k).t;
        dataStruct(k).z_r = NaN(size(dataStruct(k).z));
        dataStruct(k).d_r = NaN(size(dataStruct(k).d));
        dataStruct(k).t_r_smooth = dataStruct(k).t_r;
        dataStruct(k).z_r_smooth = NaN(size(dataStruct(k).z));
        dataStruct(k).d_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).tip_rep_pos = 1;
        dataStruct(k).tip_rep_pos_smooth = 1;
        tip_rep_pos_all(k) = 1;
        dSmoothMinAll(k) = 1;
        dataStruct(k).F_r = NaN(size(dataStruct(k).d));
        dataStruct(k).h_r = NaN(size(dataStruct(k).d));
        dataStruct(k).F_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).h_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).z_max_ind = length(dataStruct(k).z);
        dataStruct(k).z_max_ind_smooth = length(dataStruct(k).z);
        dataStruct(k).dSmoothMin = 1;
        dataStruct(k).F_r_log = NaN(size(dataStruct(k).d));
        dataStruct(k).t_r_log = NaN(size(dataStruct(k).d));
        
%         if createAverage
%             n_rows = size(dataStruct,2);
%             v_approach_temp = round(v_approach,2,'significant');
%             v_unique_temp = (unique(v_approach_temp));
%             for kk = 1:length(v_unique_temp)
%                 dataStruct(n_rows+kk).z_average = NaN(size(dataStruct(k).z));
%                 dataStruct(n_rows+kk).d_average = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).t_average = dataStruct(k).t;
%                 dataStruct(n_rows+kk).t_r = NaN(size(dataStruct(k).t));
%                 dataStruct(n_rows+kk).z_r = NaN(size(dataStruct(k).z));
%                 dataStruct(n_rows+kk).d_r = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).t_r_smooth = dataStruct(k).t;
%                 dataStruct(n_rows+kk).z_r_smooth = NaN(size(dataStruct(k).z));
%                 dataStruct(n_rows+kk).d_r_smooth = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).tip_rep_pos = 1;
%                 dataStruct(n_rows+kk).tip_rep_pos_smooth = 1;
%                 tip_rep_pos_all(n_rows+kk) = 1;
%                 dSmoothMinAll(n_rows+kk) = 1;
%                 dataStruct(n_rows+kk).F_r = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).h_r = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).F_r_smooth = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).h_r_smooth = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).z_max_ind = length(dataStruct(k).z);
%                 dataStruct(n_rows+kk).z_max_ind_smooth = length(dataStruct(k).z);
%                 dataStruct(n_rows+kk).dSmoothMin = 1;
%                 dataStruct(n_rows+kk).F_r_log = NaN(size(dataStruct(k).d));
%                 dataStruct(n_rows+kk).t_r_log = NaN(size(dataStruct(k).d));
%             end
%         end
        
        continue;

    end
    
    dt = dataStruct(k).dt;    
    
    % Create a "clean" repulsive time array
    t_rep = linspace(0,(n_offset-1)*dt,n_offset);
    t_rep_smooth = linspace(0,(n_offset_smooth-1)*dt,n_offset_smooth);

    % Store the repulsive z-sensor and deflection for normal and smooth
    % cases.
    z_rep = dataStruct(k).z_corrected(dSmoothMin:end);
    d_rep = dataStruct(k).d_corrected(dSmoothMin:end);
    z_rep_smooth = dataStruct(k).z_smooth(dSmoothMin:end);
    d_rep_smooth = dataStruct(k).d_smooth(dSmoothMin:end);
    
    % Store the tip position
    tip_rep = d_rep;
    tip_rep_smooth = d_rep_smooth;
    
    % Search for the repulsive portion of the data starting from either the
    % first index and moving to the end ('forward') or starting at the end
    % and walking backward ('reverse') and looking for where the data is
    % appreciably above the approach data. Alternately, perform a Legendre
    % transformation on the data with a moving window (modified Wavelet
    % Transform) to remove noise and make the repulsive regime obvious.
    if strcmp(findRep,'legendre')
    
        % Use the Legendre-Wavelet approach
        data = d_rep.*k_cantilever(k);
        if ~includeRetract
            % Make the data symmetric
            data = horzcat(data, flip(data));
        end
        
        n_data = numel(data);
        n_legs = 10;
        n_window = 10;

        data_split = NaN((n_data-n_window),n_window);
        ai = 0;
        bi = n_window-1;
        for ii = 1:size(data_split,1)
            ai = ai+1;
            bi = bi+1;
            if bi < numel(data)
                data_split(ii,:) = data(ai:bi);
            else
                break;
            end
        end

        L_out = NaN(n_window,n_legs);
        x_temp = linspace(-1,1,size(data_split,2));
        L_func = legendre(n_legs-1,x_temp);
        L_smooth = NaN(size(L_out,1)-1,1);

        for ii = 1:size(data_split,1)
            for jj = 1:n_legs
                % Legendre Transform
                L_out(ii,jj) = dot(data_split(ii,:),L_func(jj,:));
            end
            if ii > 1
                L_smooth(ii-1) = dot(L_out(ii,:),L_out(ii-1,:));
            end
        end

        data = (data(1:numel(data)/2));
        L_smooth = (L_smooth(1:ceil(numel(L_smooth)/2)));

        temp = 5*max(abs(L_smooth(1:round(numel(L_smooth)/20))));
        upperLim = mean(L_smooth(1:round(numel(L_smooth)/20)))+temp;
        idx = find(L_smooth>upperLim,1,'first')-1;

        if ~isempty(idx)
            tip_rep_pos = idx;
            if idx-delay > 0
                tip_rep_pos_smooth = idx-delay;
            else
                tip_rep_pos_smooth = 2;
            end
        else
            tip_rep_pos = 2;
            tip_rep_pos_smooth = 2;
        end

%         % Visualize the contact point estimate!
%         try 
%             figure(tfig);
%             clf
%             tiledlayout(2,1);
%         catch
%             tfig = figure;
%             tiledlayout(2,1);
%         end
% 
%         nexttile
%         plot(1:numel(L_smooth),L_smooth,'r-','linewidth',3)
%         hold on
%         scatter(idx,L_smooth(idx),'bo','linewidth',3)
%         plot(1:numel(L_smooth),upperLim*ones(size(1:numel(L_smooth))),'-b','linewidth',3)
%         xlabel('Bin Number')
%         ylabel('Magnitude')
%         xlim([1 numel(L_smooth)])
%         title('Smooth Bin Data')
%         hold off
% 
%         nexttile
%         plot(1:numel(data),data,'r-','linewidth',3)
%         hold on
%         scatter(idx,data(idx),'bo','linewidth',3)
%         xlabel('Data Index')
%         ylabel('Force [N]')
%         xlim([1 numel(data)])
%         title('Original Force Data')
%         hold off
% 
%         disp('pause');
        
    else
    
        % Legacy Method
        if strcmp(findRep,'forward')
            tip_rep_pos = find(tip_rep>0,1);                                   % Find first position above 0
            if isempty(tip_rep_pos) || tip_rep_pos < 2
                tip_rep_pos = 2;
            end
            
            tip_rep_pos_smooth = find(tip_rep_smooth>0,1);                     % Find first position above 0
            if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth < 2
                tip_rep_pos_smooth = 2;
            end
        elseif strcmp(findRep,'reverse')
            tip_rep_pos = (length(tip_rep) - find(flip(tip_rep)<0,1));       % Find last position above 0
            if isempty(tip_rep_pos) || tip_rep_pos < 2
                tip_rep_pos = 2;
            end
            
            tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flip(tip_rep_smooth)<0,1));   % Find last position above 0
            if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth < 2
                tip_rep_pos_smooth = 2;
            end
        end
        
    end
    
    % Store the repulsive force application portion of the dataset.
    dataStruct(k).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos-1);
    dataStruct(k).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos-1);
    dataStruct(k).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos-1);
    dataStruct(k).t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth-1);
    dataStruct(k).z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth-1);
    dataStruct(k).d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth-1);
    
    dataStruct(k).tip_rep_pos = tip_rep_pos;
    dataStruct(k).tip_rep_pos_smooth = tip_rep_pos_smooth;
    
    % Store for benchmarking later. This will help us decide which region
    % of the curves we can actually "average" together, since we want to
    % use data where all files that are averaged have some contribution
    % (i.e. where they overlap).
    tip_rep_pos_all(k) = tip_rep_pos;
    dSmoothMinAll(k) = dSmoothMin;
    
    % Calculate Force and Indentation during Repulsive Portion
    dataStruct(k).F_r = dataStruct(k).d_r .* k_cantilever(k); % Calculate Force
    dataStruct(k).h_r = (dataStruct(k).z_r - dataStruct(k).d_r); % Calculate Indentation
    
    dataStruct(k).F_r_smooth = dataStruct(k).d_r_smooth .* k_cantilever(k); % Calculate Smooth Force
    dataStruct(k).h_r_smooth = dataStruct(k).z_r_smooth - dataStruct(k).d_r_smooth; % Calculate Smooth Indentation
    
    % Store our critical indices
    dataStruct(k).z_max_ind = z_max_ind;
    dataStruct(k).z_max_ind_smooth = z_max_ind-delay;
    dataStruct(k).dSmoothMin = dSmoothMin;

    % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
    tr = dataStruct(k).dt;
    st = dataStruct(k).t(z_max_ind);
    
    F_log = [];
    t_log = [];
    
    F_log = log_scale(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
    t_log = log_scale(dataStruct(k).t_r,dataStruct(k).t_r,tr,st);
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(k).F_r_log = F_log;
    dataStruct(k).t_r_log = t_log;

end

clearvars k point_number_temp run_number_temp

% Determine how many load levels we are considering
v_approach = round(v_approach,2,'significant');
v_unique = (unique(v_approach));
r_tip_array = zeros(size(v_unique));
nu_sample_array = zeros(size(v_unique));

n_rows = size(dataStruct,2);

% Pre-Processing for Averaged Data of all load levels.
if createAverage

    if (size(dataStruct,2) > 1)
        N = n_rows;

        for k = 1:n_rows
            [minVal(k),~] = min(dataStruct(k).t_approach);
            [maxVal(k),~] = max(dataStruct(k).t_approach);
            dtVal(k) = dataStruct(k).dt;
        end

        if length(v_unique) == 1

            % Find files corresponding to current velocity
            [startNum,~] = max(minVal);
            [minRepInd,minRepFile] = min(tip_rep_pos_all + dSmoothMinAll);
            [endNum,~] = min(maxVal);

            xi = (startNum:median(dtVal):endNum)';  % Create Vector Of Common time values
            di = [];
            zi = [];
            ti = [];
            dataLengths = [];
            currentFile = [];
            forceLimits = [];
            timeLimits = [];

            r_tip_array(1) = mode(cell2mat({dataStruct(1:n_rows).r_tip}));
            nu_sample_array(1) = mode(cell2mat({dataStruct(1:n_rows).nu_sample}));

            % When only one velocity is present
            for k = 1:n_rows
                if k == minRepFile
                    adjustInd = 1;
                else
                    adjustInd = 0;
                end
                shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin-1) - ...
                                dataStruct(minRepFile).t_approach(minRepInd-adjustInd));

                [tempz, ~] = unique(dataStruct(k).t_approach - shiftVal);

                di(:,k) = interp1(tempz, dataStruct(k).d_corrected,...
                    xi(:), 'linear', NaN); % Interploate deflection to new �x� Values

                zi(:,k) = interp1(tempz, dataStruct(k).z_corrected,...
                    xi(:), 'linear', NaN); % Interploate z-sensor to new �x� Values

                % Create an associated time array for averaging
                ti(:,k) =  xi;

                % Hold on to the data lengths so we can keep track
                % of which files are worst and should be ignored.
                if isempty(find(isnan(zi(k,:)),1,'first'))
                    dataLengths(k) = size(zi(k,:),2);
                    forceLimits(k) = max(di(k,:));
                    timeLimits(k) = max(tempz);
                else
                    dataLengths(k) = find(isnan(zi(k,:)),1,'first');
                    forceLimits(k) = max(di(k,:));
                    timeLimits(k) = max(tempz);
                end
                currentFile(k) = k;

            end

            % Interpolate according to the variance in the time array that
            % comes from averaging the results.
            t_interp = (1:size(ti,2))*median(dtVal);
            z_interp = interp1(mean(ti,2), mean(zi,2),...
                    t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values
            d_interp = interp1(mean(ti,2), mean(di,2),...
                    t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values

            nanCheck = isnan(reshape(t_interp,1,length(t_interp))) +...
                isnan(reshape(z_interp,1,length(z_interp))) +...
                isnan(reshape(d_interp,1,length(d_interp)));
            nanCheck(nanCheck ~= 0) = 1;

            t_interp(nanCheck == 1) = [];
            z_interp(nanCheck == 1) = [];
            d_interp(nanCheck == 1) = [];

            if size(z_interp,2) < size(z_interp,1)
                dataStruct(n_rows+1).z_average = z_interp';
            else
                dataStruct(n_rows+1).z_average = z_interp;
            end

            if size(d_interp,2) < size(d_interp,1)
                dataStruct(n_rows+1).d_average = d_interp';
            else
                dataStruct(n_rows+1).d_average = d_interp;
            end

            if size(t_interp,2) < size(t_interp,1)
                dataStruct(n_rows+1).t_average = t_interp';
            else
                dataStruct(n_rows+1).t_average = t_interp;
            end

            % Sanity check --- should be the same as the datasets' dt
            dataStruct(n_rows+1).dt = mean(diff(dataStruct(n_rows+1).t_average));

        else

            % Loop through unique velocities
            for j = 1:length(v_unique)

                % Find files corresponding to current velocity
                velInd = zeros(1,N);
                velInd(v_approach == v_unique(j)) = 1;

                [startNum,~] = max(minVal(velInd==1));
                [minRepInd,temp] = min(tip_rep_pos_all(velInd==1) + dSmoothMinAll(velInd==1));
                fileInds = find(velInd);
                minRepFile = fileInds(temp);
                [endNum,~] = min(maxVal(velInd==1));

                xi = (startNum:median(dtVal):endNum)';  % Create Vector Of Common time values
                di = [];
                zi = [];
                ti = [];
                dataLengths = [];
                currentFile = [];
                forceLimits = [];
                timeLimits = [];

                r_tip_array(j) = mode(cell2mat({dataStruct(v_approach(:) == v_unique(j)).r_tip}));
                nu_sample_array(j) = mode(cell2mat({dataStruct(v_approach(:) == v_unique(j)).nu_sample}));

                if sum(v_approach(:) == v_unique(j)) > 1
                    % Loop through all files
                    for k = 1:n_rows
                        % Grab the tip radius for this velocity and 

                        % Check if this file is relevant for this specific
                        % averaging operation
                        if v_approach(k) == v_unique(j)
                            if k == minRepFile
                                adjustInd = 1;
                            else
                                adjustInd = 0;
                            end
                            shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin-1) - ...
                                dataStruct(minRepFile).t_approach(minRepInd-adjustInd));

                            [tempz, ~] = unique(dataStruct(k).t_approach - shiftVal);

                            di = horzcat(di, interp1(tempz, dataStruct(k).d_corrected,...
                                xi(:), 'linear', NaN)); % Interploate deflection to new �x� Values

                            zi = horzcat(zi, interp1(tempz, dataStruct(k).z_corrected,...
                                xi(:), 'linear', NaN)); % Interploate z-sensor to new �x� Values

                            % Create an associated time array for averaging
                            ti =  horzcat(ti, xi);

                            % Hold on to the data lengths so we can keep track
                            % of which files are worst and should be ignored.
                            if isempty(find(isnan(zi(:,end)),1,'first'))
                                dataLengths = horzcat(dataLengths, size(zi(:,end),2));
                                forceLimits = horzcat(forceLimits, max(di(:,end)));
                                timeLimits = horzcat(timeLimits, max(tempz));
                            else
                                dataLengths = horzcat(dataLengths, find(isnan(zi(:,end)),1,'first'));
                                forceLimits = horzcat(forceLimits, max(di(:,end)));
                                timeLimits = horzcat(timeLimits, max(tempz));
                            end
                            currentFile = horzcat(currentFile, k);

                        else
                            continue;
                        end

                    end
                else
                    % There is only one file for this velocity, so treat it
                    % like a single file would be then skip to the next
                    % iteration.
                    velIndRef = find(velInd);

                    if size(dataStruct(velIndRef).z_corrected,2) < size(dataStruct(velIndRef).z_corrected,1)
                        dataStruct(n_rows+j).z_average = dataStruct(velIndRef).z_corrected';
                    else
                        dataStruct(n_rows+j).z_average = dataStruct(velIndRef).z_corrected;
                    end

                    if size(dataStruct(velIndRef).d_corrected,2) < size(dataStruct(velIndRef).d_corrected,1)
                        dataStruct(n_rows+j).d_average = dataStruct(velIndRef).d_corrected';
                    else
                        dataStruct(n_rows+j).d_average = dataStruct(velIndRef).d_corrected;
                    end

                    if size(dataStruct(velIndRef).t_approach,2) < size(dataStruct(velIndRef).t_approach,1)
                        dataStruct(n_rows+j).t_average = dataStruct(velIndRef).t_approach';
                    else
                        dataStruct(n_rows+j).t_average = dataStruct(velIndRef).t_approach;
                    end

                    dataStruct(n_rows+j).dt = dataStruct(velIndRef).dt;
                    continue;

                end

                % Average this load level's curves.
                % Interpolate according to the variance in the time array that
                % comes from averaging the results.
                t_interp = (1:size(ti,2))*median(dtVal);
                z_interp = interp1(mean(ti,2), mean(zi,2),...
                        t_interp, 'linear', NaN); % Interploate To New Time Values
                d_interp = interp1(mean(ti,2), mean(di,2),...
                        t_interp, 'linear', NaN); % Interploate To New Time Values

                nanCheck = isnan(reshape(t_interp,length(t_interp),1)) +...
                    isnan(reshape(z_interp,length(z_interp),1)) +...
                    isnan(reshape(d_interp,length(d_interp),1));
                nanCheck(nanCheck ~= 0) = 1;

                t_interp(nanCheck == 1) = [];
                z_interp(nanCheck == 1) = [];
                d_interp(nanCheck == 1) = [];

                if size(z_interp,2) < size(z_interp,1)
                    dataStruct(n_rows+j).z_average = z_interp';
                else
                    dataStruct(n_rows+j).z_average = z_interp;
                end

                if size(d_interp,2) < size(d_interp,1)
                    dataStruct(n_rows+j).d_average = d_interp';
                else
                    dataStruct(n_rows+j).d_average = d_interp;
                end

                if size(t_interp,2) < size(t_interp,1)
                    dataStruct(n_rows+j).t_average = t_interp';
                else
                    dataStruct(n_rows+j).t_average = t_interp;
                end

                % Sanity check --- should be the same as the datasets' dt
                dataStruct(n_rows+j).dt = mean(diff(dataStruct(n_rows+1).t_average));

            end
        end

    else
        % There is only one file being analyzed, so simply copy over the data
        % to the "average" row.
        dataStruct(n_rows+1).z_average = dataStruct(1).z_corrected;
        dataStruct(n_rows+1).d_average = dataStruct(1).d_corrected;
        dataStruct(n_rows+1).t_average = dataStruct(1).t_approach;
        dataStruct(n_rows+1).dt = dataStruct(1).dt;
        r_tip_array(1) = (cell2mat({dataStruct(1).r_tip}));
        nu_sample_array(1) = (cell2mat({dataStruct(1).nu_sample}));
    end
    
    for i = 1:length(v_unique)

        % Store the tip radius in this row for reference later
        dataStruct(n_rows+i).r_tip = r_tip_array(i);
        dataStruct(n_rows+i).nu_sample = nu_sample_array(i);

        % Find the approach portion of the data. This is only really necessary
        % if the average data happens to go past the limit of all the datasets
        % for some reason.
        [~, z_max_ind] = max(dataStruct(n_rows+i).z_average);

        if ~includeRetract
            dataStruct(n_rows+i).z_approach = dataStruct(n_rows+i).z_average(1:z_max_ind);
            dataStruct(n_rows+i).d_approach = dataStruct(n_rows+i).d_average(1:z_max_ind);
        else
            F_temp = dataStruct(n_rows+i).d_average(z_max_ind:end) .* mean(k_cantilever);
            if ~isempty(F_temp) && length(F_temp) > 1
                non_contact_ind = find(F_temp < 0,1);
                if isempty(non_contact_ind) non_contact_ind = length(F_temp)-1; end
            else
                non_contact_ind = 0;
            end
            z_max_ind = non_contact_ind+z_max_ind;
            dataStruct(n_rows+i).z_approach = dataStruct(n_rows+i).z_average(1:z_max_ind);
            dataStruct(n_rows+i).d_approach = dataStruct(n_rows+i).d_average(1:z_max_ind);
        end

        % Filter and Shift z_sensor data
        switch filterType
            case 'butter'
                % Create the butterworth
                [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(n_rows+i).dt)),'low'); % This makes a lowpass filter
                d_approach_smooth = (filter(b,a,dataStruct(n_rows+i).d_approach)); % Next, apply the filter

                [~, dSmoothMin] = min(d_approach_smooth);   
                z_approach_smooth = dataStruct(n_rows+i).z_approach;

            case {'IIR','FIR'}
                % Use an IIR or FIR filter on the data
                Fs = 1/(dataStruct(n_rows+i).dt); 
                Fstop = 100*( (round((v_unique(i)),2,'significant')...
                    /round(max(v_unique),2,'significant'))...
                    /dataStruct(n_rows+i).dt );
                if Fstop >= 1/(3*dataStruct(n_rows+i).dt)
                    Fstop = 1/(3*dataStruct(n_rows+i).dt);
                elseif Fstop < 1/(10*dataStruct(n_rows+i).dt)
                    Fstop = 1/(10*dataStruct(n_rows+i).dt);
                end
                Fpass = Fstop*0.01;
                Rp = 0.01;
                Astop = 80;
                LPF = dsp.LowpassFilter('SampleRate',Fs, ...
                                         'FilterType',filterType, ...
                                         'PassbandFrequency',Fpass, ...
                                         'StopbandFrequency',Fstop, ...
                                         'PassbandRipple',Rp, ...
                                         'StopbandAttenuation',Astop);
                delay = floor(mean(grpdelay(LPF)));

                d_approach_smooth = LPF(dataStruct(n_rows+i).d_approach);         % Smooth with IIR filter

                % Correct filter delay
                sf = d_approach_smooth;
                sf(1:delay) = [];

                [~, dSmoothMin] = min(sf);
                d_approach_smooth = sf;
                z_approach_smooth = dataStruct(n_rows+i).z_approach((delay+1):end);

            case 'none'
                % Apply no filtering
                d_approach_smooth = dataStruct(n_rows+i).d_approach;
                delay = 0;

                [~, dSmoothMin] = min(d_approach_smooth);   
                z_approach_smooth = dataStruct(n_rows+i).z_approach;

        end

        t_rep = dataStruct(n_rows+i).t_average(dSmoothMin:z_max_ind);
        z_rep = dataStruct(n_rows+i).z_average(dSmoothMin:z_max_ind);
        d_rep = dataStruct(n_rows+i).d_average(dSmoothMin:z_max_ind);
        t_rep_smooth = dataStruct(n_rows+i).t_average(dSmoothMin:(z_max_ind-delay));
        z_rep_smooth = z_approach_smooth(dSmoothMin:(z_max_ind-delay));
        d_rep_smooth = d_approach_smooth(dSmoothMin:(z_max_ind-delay));

        tip_rep = d_rep;
        tip_rep_smooth = d_rep_smooth;
        if strcmp(findRep,'forward')
            tip_rep_pos = find(tip_rep>0,1);                                   % Find first position above 0
            tip_rep_pos_smooth = find(tip_rep_smooth>0,1);                     % Find first position above 0
        elseif strcmp(findRep,'reverse')
            tip_rep_pos = (length(tip_rep) - find(flip(tip_rep)<0,1));       % Find last position above 0
            if isempty(tip_rep_pos) || tip_rep_pos == 0
                tip_rep_pos = 1;
            end

            tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flip(tip_rep_smooth)<0,1));       % Find last position above 0
            if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth == 0
                tip_rep_pos_smooth = 1;
            end
        end

        % Find the repulsive portion (force application) region
        dataStruct(n_rows+i).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos);
        dataStruct(n_rows+i).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos);
        dataStruct(n_rows+i).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos);
        t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth);
        z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth);
        d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth);

        % Calculate Force and Indentation during Repulsive Portion
        dataStruct(n_rows+i).F_r = dataStruct(n_rows+i).d_r .* mean(k_cantilever);
        dataStruct(n_rows+i).h_r = dataStruct(n_rows+i).z_r - dataStruct(n_rows+i).d_r; % Calculate Indentation

        dataStruct(n_rows+i).z_max_ind = z_max_ind;
        dataStruct(n_rows+i).z_max_ind_smooth = z_max_ind - delay;
        dataStruct(n_rows+i).dSmoothMin = dSmoothMin;

        % Create smooth Force and Indentation
        k = 1;
        k_avg = 0;
        while k_avg == 0
            % Check if this file is relevant for finding cantilever stiffness
            if v_approach(k) == v_unique(i)
                k_avg = k_cantilever(i);
                continue;
            else
                k = k+1;
                continue;
            end
        end

        dataStruct(n_rows+i).t_r_smooth = t_r_smooth;
        dataStruct(n_rows+i).z_r_smooth = z_r_smooth;
        dataStruct(n_rows+i).d_r_smooth = d_r_smooth;
        dataStruct(n_rows+i).F_r_smooth = d_r_smooth .* k_avg; % Calculate Smooth Force
        dataStruct(n_rows+i).h_r_smooth = z_r_smooth - d_r_smooth; % Calculate Smooth Indentation

        if removeNegatives
            % Original
            toRemove = (dataStruct(n_rows+i).h_r <= 0 | dataStruct(n_rows+i).F_r <= 0);
            dataStruct(n_rows+i).t_r(toRemove) = [];
            dataStruct(n_rows+i).z_r(toRemove) = [];
            dataStruct(n_rows+i).d_r(toRemove) = [];
            dataStruct(n_rows+i).F_r(toRemove) = [];
            dataStruct(n_rows+i).h_r(toRemove) = [];

            % Smooth
            toRemove = (dataStruct(n_rows+i).h_r_smooth <= 0 | dataStruct(n_rows+i).F_r_smooth <= 0);
            dataStruct(n_rows+i).t_r_smooth(toRemove) = [];
            dataStruct(n_rows+i).z_r_smooth(toRemove) = [];
            dataStruct(n_rows+i).d_r_smooth(toRemove) = [];
            dataStruct(n_rows+i).F_r_smooth(toRemove) = [];
            dataStruct(n_rows+i).h_r_smooth(toRemove) = [];
        end

        % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
        tr = dataStruct(n_rows+i).dt;
        st = dataStruct(n_rows+i).t_r(end);

        F_log = [];
        t_log = [];

        F_log = log_scale(dataStruct(n_rows+i).F_r,dataStruct(n_rows+i).t_r,tr,st);
        t_log = log_scale(dataStruct(n_rows+i).t_r,dataStruct(n_rows+i).t_r,tr,st);

        % Save the Log-Form of the tip force (F) and time array (t_r)
        dataStruct(n_rows+i).F_r_log = F_log;
        dataStruct(n_rows+i).t_r_log = t_log;
    end
    
end

end

