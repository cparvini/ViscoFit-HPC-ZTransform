function [] = globalMapClustering(originalPath,N_workers,clusterTarget,varargin)
%GLOBALMAPCLUSTERING Perform K-Medoids Clustering of Many QI Map
%   This function takes in a path argument, the target variable for
%   clustering with k-medoids, and a variable number of arguments used to
%   correct/shift the map representations. Note that the only non-boolean
%   option for varargin is evalPt, which is the frequency (in Hz) to use
%   for plotting the frequency-dependent properties from the QI map. The
%   full spectrum of frequencies available for each pixel are used in the
%   clustering process, and evalPt is exclusively for visualization
%   purposes. This function has a pair, singleMapClustering(), which will
%   perform nearly the same analysis except it only considers a single QI
%   map at once.
%
%   By searching recursively in a directory, this file creates a large
%   array of clustering data that includes all of the cell pixels from
%   every map in the subdirectories under the target. This is intended to
%   be run on QI maps for samples of the same condition (e.g. cell line and
%   stage, treatment, etc.) such that the clustering is less susceptible to
%   single-map noise or feature issues.

% User-Defined Settings
hideSubstrate = true;
fillPixels = true;
logSteps = true;
plotIndentation = true;
evalPt = 1000;
n_reps = 100; % number of clustering replicates
maxK = 10; % Max number of cluster bins
manualK = []; % Specific number of bins to use
if nargin > 1
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        hideSubstrate = varargin{i};
                    end
                case 2
                    if ~isempty(varargin{i})
                        fillPixels = varargin{i};                        
                    end
                case 3
                    if ~isempty(varargin{i})
                        logSteps = varargin{i};                        
                    end
                case 4
                    if ~isempty(varargin{i})
                        plotIndentation = varargin{i};                        
                    end
                case 5
                    if ~isempty(varargin{i})
                        evalPt = varargin{i};                        
                    end
                case 6
                    if ~isempty(varargin{i})
                        n_reps = varargin{i};                        
                    end
                case 7
                    if ~isempty(varargin{i})
                        maxK = varargin{i};                        
                    end
                case 8
                    if ~isempty(varargin{i})
                        manualK = varargin{i};                        
                    end
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

% Permanent Settings
errortype = 'sse';
figX = 0;
figY = 0;
maxCol = 5;
maxwid = get(0,'screensize');
maxwid = maxwid(3);
mapColorName = 'turbo';
mapEdgeCol = 'none';
climMax = 2e5; % Pa
climHeight = 15e-6; % meters, the JPK Nanowizard has a 15um piezo 
climInd = 1000e-9; % meters
trimHeight = 100e-9;
dFreq = 200; % Hz, step size between frames, if discrete
n_steps = 50; % frames, number of frames per order of magnitude
n_datapoints = 10;

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~iscell(originalPath)
    % Perform global clustering on one directory
    Folders = {originalPath};
else
    % The user has given multiple 
    Folders = originalPath;
end

% Define error functions
sse_global = @(data,model) sum((data-model).^2,'all');
mse_global = @(data,model,n) sum((data-model).^2,'all')./(length(data)-n);

% Clear Previous Parpool
if ~isempty(gcp('nocreate'))
   % Get the current pool
    poolobj = gcp('nocreate');
    if poolobj.NumWorkers ~= N_workers
        delete(poolobj);
        % Make a fresh pool
        if N_workers == 1 
            poolobj = [];
            parallelSet = false;
        elseif N_workers > 1
            poolobj = parpool(N_workers,'IdleTimeout', 120);
            parallelSet = true;
        else
            warning('Warning: you provided an invalid value for "N_workers". Using system default configuration.');
            poolobj = parpool('IdleTimeout', 120);
            parallelSet = true;
        end
    else
        if poolobj.NumWorkers > 1
            parallelSet = true;
        end
    end
else
    % Make a fresh pool
    if N_workers == 1 
        poolobj = [];
        parallelSet = false;
    elseif N_workers > 1
        poolobj = parpool(N_workers,'IdleTimeout', 120);
        parallelSet = true;
    else
        warning('Warning: you provided an invalid value for "N_workers". Using system default configuration.');
        poolobj = parpool('IdleTimeout', 120);
        parallelSet = true;
    end
end
    
% Begin looping through the directories or files
for i_dir = 1:length(Folders)

    % Issue wrapper which allows script to continue if there is an empty
    % directory, or other issue during processing.
    try
    
        % Search recursively (using "**") for our zTransform results files
        path = Folders{i_dir};
        Files = dir(fullfile(path, '**','*Results*zTransform*.mat'));

        if isempty(Files)
            error('The directory you selected and its subdirectories do not contain a Z-Transform QI map. Please verify your results file is in or under that directory and the filename contains "zTransform".');
        end

        fileLabels = cell(numel(Files),1);
        for j = 1:numel(Files)
            temp = strsplit(Files(j).name,{'-','_','.'});
            idx = find(contains(lower(temp),{'fitresults','mapresults'}),1);
            fileLabels{j} = strjoin(temp([1 idx-1]),'-');
        end
        
        % Since we are analyzing many maps, we must create our frequency
        % array ahead of time.
        fprintf('\nCreating a frequency vector...');

        vars = whos('-file',[Files(1).folder filesep Files(1).name]);
        varNames = {vars.name};
        mapsizes = NaN(numel(Files),1);
        for j = 1:numel(varNames)

            if ~contains(varNames{j},'zTransform')
                continue;
            end

            minFreq = Inf;
            maxFreq = 0;

            for j_dir = 1:length(Files)

                % Load the mat file
                resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
                freqMapTemp = resultsStruct.(varNames{j}).frequencyMap;
                timesCellTemp = resultsStruct.(varNames{j}).ViscoClass.times_cell;
                clearvars resultsStruct
                
                mapsizes(j_dir) = numel(freqMapTemp);

                for k_pixels = 1:numel(freqMapTemp)
                    tempf = freqMapTemp{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = min(tempf,[],'omitnan');
                    if any([isnan(temp),isempty(temp)]) || (numel(timesCellTemp{k_pixels})<n_datapoints)
                        continue;
                    end
                    if temp < minFreq
                        minFreq = temp;
                    end
                    tempf = freqMapTemp{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = max(tempf,[],'omitnan');
                    if isnan(temp)|| isempty(temp)
                        continue;
                    end
                    if temp > maxFreq
                        maxFreq = temp;
                    end
                end
                
            end

            % Count orders of 10
            temp = minFreq;
            tempf = minFreq;
            tempmax = 10^ceil(log10(minFreq));
            magList = [];
            if ~logSteps
                while temp < maxFreq
                    tempf = 10.^( ( log10(temp) ) );
                    magList = horzcat(magList,tempf);
                    temp = temp + dFreq;
                end
            else
                while temp < maxFreq
                    if temp == minFreq
                        dFreq = 10^(ceil(log10(temp)))/n_steps; 
                    end

                    if temp >= tempmax
                        tempmax = temp*10;
                        dFreq = 10^(ceil(log10(temp)))/n_steps; 
                    end

                    tempf = 10.^( ( log10(temp) ) );
                    magList = horzcat(magList,tempf);
                    temp = temp + dFreq;
                end
            end
            magList = unique(magList);
            freqList = flip(magList);
            timeList = 1./freqList;

        end
        
        fprintf('\nPre-allocating arrays...');
        
        % Prepare our index records
        ai = NaN(numel(Files),1);
        bi = NaN(numel(Files),1);
        for j_dir = 1:length(Files)
            if j_dir == 1
                ai(j_dir) = 1;
            else
                ai(j_dir) = bi(j_dir-1) + 1;
            end
            bi(j_dir) = ai(j_dir) + mapsizes(j_dir) - 1;
        end
        clusteringDataGlobal = NaN(bi(end),numel(magList));
        pixelLogGlobal = NaN(bi(end),2);
        mapIDglobal = NaN(bi(end),1);
        
        fprintf('Complete!\n');
        
        % Now, loop through the files and create our dataset
        fprintf('\nPopulating our global clustering dataset...');
        for j_dir = 1:length(Files)
            
            resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
            varNames = fields(resultsStruct);

            for j = 1:numel(varNames)
                
                if ~contains(varNames{j},'zTransform')
                    continue;
                end

                if isfield(resultsStruct.(varNames{j}),'mapSize')
                    mapSize = resultsStruct.(varNames{j}).mapSize;
                else
                    mapSize = [128 128];
                end
                
                % Grab some relevant settings
                correctTilt = resultsStruct.(varNames{j}).correctTilt;
                zeroSubstrate = resultsStruct.(varNames{j}).zeroSubstrate;
                optimizeFlattening = resultsStruct.(varNames{j}).optimizeFlattening;
                
                if isfield(resultsStruct.(varNames{j}),'scanSize')
                    % We have the absolute map size!
                    scanSize = resultsStruct.(varNames{j}).scanSize;
                    xdataAbs = 0:(scanSize(1)/(mapSize(1)-1)):scanSize(1);
                    ydataAbs = flip(0:(scanSize(2)/(mapSize(2)-1)):scanSize(2));
                    [XA, YA] = meshgrid(xdataAbs,ydataAbs);
                end

                pixelHeight_cell = resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell;

                % axes meshgrid for scattering data
                xdata = 1:mapSize(1);
                ydata = flip(1:mapSize(2));
                [X, Y] = meshgrid(xdata,ydata);

                if isfield(resultsStruct.(varNames{j}),'bestParams')
                    plotModel = true;
                else
                    plotModel = false;
                end

                pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
                if correctTilt
                    temp = fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate,[],optimizeFlattening);
                    pixelHeightArray = cell2mat(temp);
                else
                    pixelHeightArray = cell2mat(pixelHeight_cell);
                end

                [minHeight,~] = min(pixelHeightArray(pixelHeightArray>0));
                substrateCutoff = minHeight + trimHeight;
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel(pixelHeightArray);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list

                % Position for the map
                xc = 1;
                yc = 0;

                if ~isfield(resultsStruct.(varNames{j}),'indMap')
                    % Do some prep to shorten the computation time and
                    % limit the number of calls to "zTransformCurve",
                    % which is slow.
                    F_hz_all = cell(size([pixelHeight_cell{:}]));
                    h_hz_all = cell(size([pixelHeight_cell{:}]));

                    for i_z = 1:numel(F_hz_all)

                        if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                            F_hz_all{i_z} = NaN;
                            h_hz_all{i_z} = NaN;
                            continue;
                        end

                        dataIn = {resultsStruct.(varNames{j}).ViscoClass.times_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.dts_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.forces_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.indentations_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.tipSize_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.nu_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.tipGeom,...
                            resultsStruct.(varNames{j}).ViscoClass.minTimescale,...
                            resultsStruct.(varNames{j}).ViscoClass.thinSample,...
                            resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell{i_z}};

                        if any(cellfun(@isempty,dataIn(1:6))) || any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                            F_hz_all{i_z} = NaN;
                            h_hz_all{i_z} = NaN;
                            continue;
                        end

                        [~,~,F_hz_all{i_z},~,h_hz_all{i_z},~,~] = zTransformCurve(dataIn,'none',0.05,resultsStruct.(varNames{j}).ViscoClass.thinSample);

                    end

                end
                
                clusteringData = NaN(numel(pixelHeightArray),numel(magList));
                pixelLog = NaN(numel(pixelHeightArray),2);

                for k_pixels = 1:numel(pixelHeightArray)

                    % Get the current pixel position
                    if xc > mapSize(1)
                        xc = 1;
                        yc = yc + 1;
                    end
                    idx_pixel = sub2ind(flip(mapSize),mapSize(2)-yc,xc);
                    pixelLog(k_pixels,:) = [mapSize(2)-yc,xc];
                    
                    if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                        xc = xc + 1;
                        continue;
                    end

                    if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                        xc = xc + 1;
                        continue;
                    end

                    if ~isfield(resultsStruct.(varNames{j}),'indMap')
                        dataIn = {resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.dts_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.tipSize_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.nu_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.tipGeom,...
                            resultsStruct.(varNames{j}).ViscoClass.minTimescale,...
                            resultsStruct.(varNames{j}).ViscoClass.thinSample,...
                            resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell{k_pixels}};

                        [~,~,F_hz,~,h_hz,~,~] = zTransformCurve(dataIn,'none',0.05,resultsStruct.(varNames{j}).ViscoClass.thinSample);
                        F_hz = abs(F_hz);
                        h_hz = abs(h_hz);
                    else
                        F_hz = abs(resultsStruct.(varNames{j}).forceMap{k_pixels});
                        h_hz = abs(resultsStruct.(varNames{j}).indMap{k_pixels});
                    end

                    % Load and perform peak correction
                    freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    [~,maxid] = max(F_hz);
                    freqAdj = freq(maxid);
                    freq = freq - freqAdj;

                    % Resample to known array of frequencies
                    ids = ((freq >= min(magList)) & (freq <= max(magList)));

                    if sum(ids,'all') < 2
                        xc = xc + 1;
                        continue;
                    end

                    if plotModel

                        harmonicSettings = struct;
                        harmonicSettings.elasticSetting = resultsStruct.(varNames{j}).elasticSetting;
                        harmonicSettings.fluidSetting = resultsStruct.(varNames{j}).fluidSetting;
                        harmonicSettings.model = resultsStruct.(varNames{j}).model;

                        % Create a frequency array
                        visco = resultsStruct.(varNames{j}).ViscoClass;
                        elasticSetting = resultsStruct.(varNames{j}).elasticSetting;
                        fluidSetting = resultsStruct.(varNames{j}).fluidSetting;

                        % Generate a frequency array in log scale
                        omega = 2.*pi.*freq;

                        % Find the best number of terms for this pixel
                        paramErrors = Inf(numel(resultsStruct.(varNames{j}).bestParams),1);
                        for k = 1:numel(resultsStruct.(varNames{j}).bestParams)

                            % Calculate the error for those parameters
                            switch lower(resultsStruct.(varNames{j}).model)
                                case 'maxwell'
                                    % Grab the parameters
                                    tempParams = resultsStruct.(varNames{j}).bestParams{k}{k_pixels};
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.forces_cell{k_pixels},...
                                                LR_Maxwell(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.forces_cell{k_pixels},...
                                                LR_Maxwell(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(tempParams)*length(visco.tipSize_cell{k_pixels}));
                                    end

                                case 'voigt'
                                    % Grab the parameters
                                    tempParams = resultsStruct.(varNames{j}).bestParams{k}{k_pixels};
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.indentations_cell{k_pixels},...
                                                LR_Voigt(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.indentations_cell{k_pixels},...
                                                LR_Voigt(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(tempParams)*length(visco.tipSize_cell{k_pixels}));
                                    end

                                case 'plr'
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.indentations_cell{k_pixels},...
                                                LR_PLR(resultsStruct.(varNames{j}).bestParams{2}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.indentations_cell{k_pixels},...
                                                LR_PLR(resultsStruct.(varNames{j}).bestParams{1}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(resultsStruct.(varNames{j}).bestParams{1}{k_pixels})*length(visco.tipSize_cell{k_pixels}));
                                    end

                                otherwise
                                    error('The model in your results structure was not recognized.')
                            end

                        end

                        % Determine the best number of arms for this pixel
                        [~,idx] = min(paramErrors);
                        bestidx = idx;
                        for k = 1:length(paramErrors)
                            if k == idx
                                continue;
                            end
                            if 100*(paramErrors(idx)-paramErrors(k))/paramErrors(idx) < 1
                                bestidx = k;                                
                            end
                        end

                        harmonicSettings.bestParams = resultsStruct.(varNames{j}).bestParams{bestidx}{k_pixels};

                        switch lower(resultsStruct.(varNames{j}).model)
                            case 'maxwell'
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_Maxwell(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.forces_cell{k_pixels},...
                                            LR_Maxwell(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.forces_cell{k_pixels},...
                                            LR_Maxwell(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(harmonicSettings.bestParams)*length(visco.tipSize_cell{k_pixels}));
                                end

                            case 'voigt'
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_Voigt(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.indentations_cell{k_pixels},...
                                            LR_Voigt(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.indentations_cell{k_pixels},...
                                            LR_Voigt(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(harmonicSettings.bestParams)*length(visco.tipSize_cell{k_pixels}));
                                end

                            case 'plr'
                                harmonicSettings.dt = dt;
                                harmonicSettings.nu_sample = mode(cell2mat(cellfun(@(x)mode(round(x,4,'significant')),nu,'UniformOutput',false)));
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_PLR(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.indentations_cell{k_pixels},...
                                            LR_PLR(resultsStruct.(varNames{j}).bestParams{2}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.indentations_cell{k_pixels},...
                                            LR_PLR(resultsStruct.(varNames{j}).bestParams{1}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(resultsStruct.(varNames{j}).bestParams{1}{k_pixels})*length(visco.tipSize_cell{k_pixels}));
                                end

                            otherwise
                                error('The model in your results structure was not recognized.')
                        end

                        modelRelaxance = modelStorage + 1j*modelLoss;

                        % Has to be from the time domain
                        % (normalization issues)
                        t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                        h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                        F_t = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};

                    else

                        % We are just plotting the data captured
                        % DIRECTLY from the z-transform method.
                        modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                        modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                        modelAngle = atand(modelLoss./modelStorage);
                        modelRelaxance = abs(resultsStruct.(varNames{j}).relaxanceMap{k_pixels});

                        % Has to be from the time domain
                        % (normalization issues)
                        t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                        h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                        F_t = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
                            
                    end
                    
                    clusterInterp = [];

                    switch clusterTarget
                        case 'force'
                            xinterp = t_t;
                            obsList = timeList;
                            clusterInterp = F_t;
                        case 'indentation'
                            xinterp = t_t;
                            obsList = timeList;
                            clusterInterp = h_t;
                        case 'storage'
                            xinterp = freq;
                            obsList = magList;
                            clusterInterp = modelStorage;
                        case 'loss'
                            xinterp = freq;
                            obsList = magList;
                            clusterInterp = modelLoss;
                        case 'angle'
                            xinterp = freq;
                            obsList = magList;
                            clusterInterp = modelAngle;
                        case 'relaxance'
                            xinterp = freq;
                            obsList = magList;
                            clusterInterp = modelRelaxance;
                    end

                    try
                        % Resample to known array of frequencies
                        obsOut = interp1(xinterp,clusterInterp,obsList,'makima',...
                            NaN);
                        
                        % SMOOTHING COULD GO HERE!
                        
                        clusteringData(k_pixels,:) = obsOut;
                    catch
                        % Do nothing
                    end
                    
                    xc = xc + 1;

                end
                
                if fillPixels
                    % Perform Interpolation of non-viable pixels
                    % (nothing for now, design something to average the
                    % time series if possible
                    
                end
                
                % Concatenate our Clustering Data from this map and save
                % the indices for this particular file
                pixelLogGlobal(ai(j_dir):bi(j_dir),:) = pixelLog;
                mapIDglobal(ai(j_dir):bi(j_dir)) = j_dir;
                clusteringDataGlobal(ai(j_dir):bi(j_dir),:) = clusteringData;

            end
            
            clearvars resultsStruct mapSize

        end
        
        fprintf('Complete!\n');
        
        fprintf('\nClustering the global dataset...');

        % Perform some data cleaning. Columns where we don't have enough
        % observations compared to bins are a PROBLEM. We have to remove
        % these before analyzing.
        dataVerify = ~isnan(clusteringDataGlobal);
        goodCols = sum(dataVerify,1);
        
        if isempty(manualK)
            % Search for the optimal number of bins
            idRem = goodCols < maxK;

            magList(idRem) = [];
            freqList(idRem) = [];
            timeList(idRem) = [];
            clusteringDataGlobal(:,idRem) = [];

            % Perform the clustering
            opts = statset('UseParallel',parallelSet,...
                'MaxIter',1e2,...
                'Display','off');
            tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
                'Distance',@dtwf,...
                'Replicates',n_reps,...
                'OnlinePhase','on',...
                'Algorithm','large');
            ids = all(isnan(clusteringDataGlobal),2); % Find excluded pixels
            clusteringDataGlobal(ids,:) = [];
            pixelLogGlobal(ids,:) = [];
            mapIDglobal(ids) = [];

            % Try all of the cluster configurations. We start with 3
            % clusters to account for any substrate pixels that were
            % not trimmed. This is not ideal, but unavoidable with
            % automatic large-scale analysis.
            eva = evalclusters(clusteringDataGlobal,tempfunc,'CalinskiHarabasz',...
                'klist',((3-hideSubstrate):maxK));

            fprintf('\nOptimal Number of Bins: %d\n\n',eva.OptimalK);

            [idxK, centroidLocs] = tempfunc(clusteringDataGlobal,eva.OptimalK);
        else
            % Use the manually specified number of bins
            idRem = goodCols < manualK;

            magList(idRem) = [];
            freqList(idRem) = [];
            timeList(idRem) = [];
            clusteringDataGlobal(:,idRem) = [];

            % Perform the clustering
            opts = statset('UseParallel',parallelSet,...
                'MaxIter',1e2,...
                'Display','off');
            tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
                'Distance',@dtwf,...
                'Replicates',n_reps,...
                'OnlinePhase','on',...
                'Algorithm','large');
            ids = all(isnan(clusteringDataGlobal),2); % Find excluded pixels
            clusteringDataGlobal(ids,:) = [];
            pixelLogGlobal(ids,:) = [];
            mapIDglobal(ids) = [];

            fprintf('\nManual Number of Bins: %d\n\n',manualK);

            [idxK, centroidLocs] = tempfunc(clusteringDataGlobal,manualK);
        end
        
        fprintf('Complete!\n');
        
        % Now, go through and save the results to each output file and make
        % our plots!
        fprintf('\nGenerating our output plots/files...');
        for j_dir = 1:length(Files)
            
            resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
            varNames = fields(resultsStruct);

            for j = 1:numel(varNames)
                
                if ~contains(varNames{j},'zTransform')
                    continue;
                end

                if isfield(resultsStruct.(varNames{j}),'mapSize')
                    mapSize = resultsStruct.(varNames{j}).mapSize;
                else
                    mapSize = [128 128];
                end
                
                % Grab some relevant settings
                correctTilt = resultsStruct.(varNames{j}).correctTilt;
                zeroSubstrate = resultsStruct.(varNames{j}).zeroSubstrate;
                optimizeFlattening = resultsStruct.(varNames{j}).optimizeFlattening;
                
                if isfield(resultsStruct.(varNames{j}),'scanSize')
                    % We have the absolute map size!
                    scanSize = resultsStruct.(varNames{j}).scanSize;
                    xdataAbs = 0:(scanSize(1)/(mapSize(1)-1)):scanSize(1);
                    ydataAbs = flip(0:(scanSize(2)/(mapSize(2)-1)):scanSize(2));
                    [XA, YA] = meshgrid(xdataAbs,ydataAbs);
                end

                pixelHeight_cell = resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell;

                % axes meshgrid for scattering data
                xdata = 1:mapSize(1);
                ydata = flip(1:mapSize(2));
                [X, Y] = meshgrid(xdata,ydata);

                if isfield(resultsStruct.(varNames{j}),'bestParams')
                    plotModel = true;
                else
                    plotModel = false;
                end

                n_plots = 3+plotIndentation;
                n_rows = ceil(n_plots/maxCol);
                n_cols = min([n_plots maxCol]);              
                mult = min([400 maxwid/n_cols]);
                figWid = mult*n_cols;
                figHeight = max([mult*n_rows figWid/n_plots]);
                
                % Clear old figures if they exist
                if ~exist('mapPlotWindow','var')
                    mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
                else
                    try
                        figure(mapPlotWindow)
                        clf
                    catch
                        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
                    end
                end

                if plotModel
                    mapType = '-Model';
                else
                    mapType = '-ZTrans';
                end

                switch clusterTarget
                    case 'force'
                        saveLabel = '-Force';
                    case 'indentation'
                        saveLabel = '-Ind';
                    case 'storage'
                        saveLabel = '-StorageMod';
                    case 'loss'
                        saveLabel = '-LossMod';
                    case 'angle'
                        saveLabel = '-LossAng';
                    case 'relaxance'
                        saveLabel = '-Relaxance';
                end
                
                plotFile = [Files(j_dir).folder filesep fileLabels{j_dir} saveLabel 'ClusteringGlobal-' varNames{j}...
                    mapType];
                
                pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
                if correctTilt
                    temp = fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate,[],optimizeFlattening);
                    pixelHeightArray = cell2mat(temp);
                else
                    pixelHeightArray = cell2mat(pixelHeight_cell);
                end

                [minHeight,~] = min(pixelHeightArray(pixelHeightArray>0));
                substrateCutoff = minHeight + trimHeight;
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel(pixelHeightArray);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list
                
%                 % Add frequency annotation
%                 pos = [figWid-175 10 150 25];
%                 str = [num2str(round(evalPt)) ' Hz'];
%                 annotation('textbox',...
%                     'Units','pixels',...
%                     'Position',pos,...
%                     'String',str,...
%                     'FitBoxToText','on',...
%                     'BackgroundColor','white',...
%                     'HorizontalAlignment','center',...
%                     'VerticalAlignment','middle',...
%                     'FontSize',16);

                % Make blank map data
                mapDataStorage = NaN(flip(mapSize));
                mapDataLoss = NaN(flip(mapSize));
                mapDataAngle = NaN(flip(mapSize));
                mapDataRelaxance = NaN(flip(mapSize));
                mapDataError = NaN(flip(mapSize));
                mapDataTerms = NaN(flip(mapSize));
                mapDataHeight = NaN(flip(mapSize));
                mapDataInd = NaN(flip(mapSize));
                mapDataForce = NaN(flip(mapSize));
                heightImg = zeros(flip(mapSize));

                % Position for the map
                xc = 1;
                yc = 0;

                if ~isfield(resultsStruct.(varNames{j}),'indMap')
                    % Do some prep to shorten the computation time and
                    % limit the number of calls to "zTransformCurve",
                    % which is slow.
                    F_hz_all = cell(size([pixelHeight_cell{:}]));
                    h_hz_all = cell(size([pixelHeight_cell{:}]));

                    for i_z = 1:numel(F_hz_all)

                        if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                            F_hz_all{i_z} = NaN;
                            h_hz_all{i_z} = NaN;
                            continue;
                        end

                        dataIn = {resultsStruct.(varNames{j}).ViscoClass.times_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.dts_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.forces_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.indentations_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.tipSize_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.nu_cell{i_z},...
                            resultsStruct.(varNames{j}).ViscoClass.tipGeom,...
                            resultsStruct.(varNames{j}).ViscoClass.minTimescale,...
                            resultsStruct.(varNames{j}).ViscoClass.thinSample,...
                            resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell{i_z}};

                        if any(cellfun(@isempty,dataIn(1:6))) || any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                            F_hz_all{i_z} = NaN;
                            h_hz_all{i_z} = NaN;
                            continue;
                        end

                        [~,~,F_hz_all{i_z},~,h_hz_all{i_z},~,~] = zTransformCurve(dataIn,'none',0.05,resultsStruct.(varNames{j}).ViscoClass.thinSample);

                    end

                end
                
                heightImg = reshape(pixelHeightArray,mapSize);

                for k_pixels = 1:numel(heightImg)

                    % Get the current pixel position
                    if xc > mapSize(1)
                        xc = 1;
                        yc = yc + 1;
                    end
                    idx_pixel = sub2ind(flip(mapSize),mapSize(2)-yc,xc);
                    
                    if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                        xc = xc + 1;
                        continue;
                    end

                    if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                        xc = xc + 1;
                        continue;
                    end

                    if ~isfield(resultsStruct.(varNames{j}),'indMap')
                        dataIn = {resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.dts_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.tipSize_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.nu_cell{k_pixels},...
                            resultsStruct.(varNames{j}).ViscoClass.tipGeom,...
                            resultsStruct.(varNames{j}).ViscoClass.minTimescale,...
                            resultsStruct.(varNames{j}).ViscoClass.thinSample,...
                            resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell{k_pixels}};

                        [~,~,F_hz,~,h_hz,~,~] = zTransformCurve(dataIn,'none',0.05,resultsStruct.(varNames{j}).ViscoClass.thinSample);
                        F_hz = abs(F_hz);
                        h_hz = abs(h_hz);
                    else
                        F_hz = abs(resultsStruct.(varNames{j}).forceMap{k_pixels});
                        h_hz = abs(resultsStruct.(varNames{j}).indMap{k_pixels});
                    end

                    % Load and perform peak correction
                    freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    [~,maxid] = max(F_hz);
                    freqAdj = freq(maxid);
                    freq = freq - freqAdj;

                    % Resample to known array of frequencies
                    ids = ((freq >= min(magList)) & (freq <= max(magList)));

                    if sum(ids,'all') < 2
                        xc = xc + 1;
                        continue;
                    end

                    if plotModel

                        harmonicSettings = struct;
                        harmonicSettings.elasticSetting = resultsStruct.(varNames{j}).elasticSetting;
                        harmonicSettings.fluidSetting = resultsStruct.(varNames{j}).fluidSetting;
                        harmonicSettings.model = resultsStruct.(varNames{j}).model;

                        % Create a frequency array
                        visco = resultsStruct.(varNames{j}).ViscoClass;
                        elasticSetting = resultsStruct.(varNames{j}).elasticSetting;
                        fluidSetting = resultsStruct.(varNames{j}).fluidSetting;

                        % Generate a frequency array in log scale
                        omega = 2.*pi.*freq;

                        % Find the best number of terms for this pixel
                        paramErrors = Inf(numel(resultsStruct.(varNames{j}).bestParams),1);
                        for k = 1:numel(resultsStruct.(varNames{j}).bestParams)

                            % Calculate the error for those parameters
                            switch lower(resultsStruct.(varNames{j}).model)
                                case 'maxwell'
                                    % Grab the parameters
                                    tempParams = resultsStruct.(varNames{j}).bestParams{k}{k_pixels};
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.forces_cell{k_pixels},...
                                                LR_Maxwell(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.forces_cell{k_pixels},...
                                                LR_Maxwell(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(tempParams)*length(visco.tipSize_cell{k_pixels}));
                                    end

                                case 'voigt'
                                    % Grab the parameters
                                    tempParams = resultsStruct.(varNames{j}).bestParams{k}{k_pixels};
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.indentations_cell{k_pixels},...
                                                LR_Voigt(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.indentations_cell{k_pixels},...
                                                LR_Voigt(tempParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(tempParams)*length(visco.tipSize_cell{k_pixels}));
                                    end

                                case 'plr'
                                    switch errortype
                                        case 'sse'
                                            paramErrors(k) = sse_global(visco.indentations_cell{k_pixels},...
                                                LR_PLR(resultsStruct.(varNames{j}).bestParams{2}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                        case 'mse'
                                            paramErrors(k) = mse_global(visco.indentations_cell{k_pixels},...
                                                LR_PLR(resultsStruct.(varNames{j}).bestParams{1}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                    numel(resultsStruct.(varNames{j}).bestParams{1}{k_pixels})*length(visco.tipSize_cell{k_pixels}));
                                    end

                                otherwise
                                    error('The model in your results structure was not recognized.')
                            end

                        end

                        % Determine the best number of arms for this pixel
                        [~,idx] = min(paramErrors);
                        bestidx = idx;
                        for k = 1:length(paramErrors)
                            if k == idx
                                continue;
                            end
                            if 100*(paramErrors(idx)-paramErrors(k))/paramErrors(idx) < 1
                                bestidx = k;                                
                            end
                        end

                        harmonicSettings.bestParams = resultsStruct.(varNames{j}).bestParams{bestidx}{k_pixels};

                        switch lower(resultsStruct.(varNames{j}).model)
                            case 'maxwell'
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_Maxwell(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.forces_cell{k_pixels},...
                                            LR_Maxwell(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.forces_cell{k_pixels},...
                                            LR_Maxwell(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.indentations_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(harmonicSettings.bestParams)*length(visco.tipSize_cell{k_pixels}));
                                end

                            case 'voigt'
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_Voigt(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.indentations_cell{k_pixels},...
                                            LR_Voigt(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.indentations_cell{k_pixels},...
                                            LR_Voigt(harmonicSettings.bestParams,visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(harmonicSettings.bestParams)*length(visco.tipSize_cell{k_pixels}));
                                end

                            case 'plr'
                                harmonicSettings.dt = dt;
                                harmonicSettings.nu_sample = mode(cell2mat(cellfun(@(x)mode(round(x,4,'significant')),nu,'UniformOutput',false)));
                                [modelStorage,modelLoss,modelAngle] = visco.harmonics_PLR(omega,harmonicSettings);
                                switch errortype
                                    case 'sse'
                                        modelErrorTime = sse_global(visco.indentations_cell{k_pixels},...
                                            LR_PLR(resultsStruct.(varNames{j}).bestParams{2}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting));
                                    case 'mse'
                                        modelErrorTime = mse_global(visco.indentations_cell{k_pixels},...
                                            LR_PLR(resultsStruct.(varNames{j}).bestParams{1}{k_pixels},visco.times_cell{k_pixels},visco.dts_cell{k_pixels},visco.forces_cell{k_pixels},visco.tipSize_cell{k_pixels},visco.nu_cell{k_pixels},visco.tipGeom,elasticSetting,fluidSetting),...
                                                numel(resultsStruct.(varNames{j}).bestParams{1}{k_pixels})*length(visco.tipSize_cell{k_pixels}));
                                end

                            otherwise
                                error('The model in your results structure was not recognized.')
                        end

                        modelRelaxance = modelStorage + 1j*modelLoss;

                        % Resample to known array of frequencies
                        mapDataStorage(idx_pixel) = interp1(freq,modelStorage,evalPt,'makima',...
                            NaN);
                        mapDataLoss(idx_pixel) = interp1(freq,modelLoss,evalPt,'makima',...
                            NaN);
                        mapDataAngle(idx_pixel) = interp1(freq,modelAngle,evalPt,'makima',...
                            NaN);
                        mapDataRelaxance(idx_pixel) = interp1(freq,modelRelaxance,evalPt,'makima',...
                            NaN);

                        % Has to be from the time domain
                        % (normalization issues)
                        t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                        h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                        F_t = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
                        mapDataInd(idx_pixel) = interp1(t_t,h_t,1/evalPt,'makima',...
                            1e-12);
                        mapDataForce(idx_pixel) = interp1(t_t,F_t,1/evalPt,'makima',...
                            1e-12);
                        
                        mapDataError(idx_pixel) = modelErrorTime;
                        mapDataTerms(idx_pixel) = bestidx;

                    else

                        % We are just plotting the data captured
                        % DIRECTLY from the z-transform method.
                        modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                        modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                        modelAngle = atand(modelLoss./modelStorage);
                        modelRelaxance = abs(resultsStruct.(varNames{j}).relaxanceMap{k_pixels});

                        % Resample to known array of frequencies                        
                        mapDataStorage(idx_pixel) = interp1(freq,modelStorage,evalPt,'makima',...
                            NaN);
                        mapDataLoss(idx_pixel) = interp1(freq,modelLoss,evalPt,'makima',...
                            NaN);
                        mapDataAngle(idx_pixel) = interp1(freq,modelAngle,evalPt,'makima',...
                            NaN);
                        mapDataRelaxance(idx_pixel) = interp1(freq,modelRelaxance,evalPt,'makima',...
                            NaN);

                        % Has to be from the time domain
                        % (normalization issues)
                        t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                        h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                        F_t = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
                        mapDataInd(idx_pixel) = interp1(t_t,h_t,1/evalPt,'makima',...
                            1e-12);
                        mapDataForce(idx_pixel) = interp1(t_t,F_t,1/evalPt,'makima',...
                            1e-12);
                            
                    end
                    
                    mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                    xc = xc + 1;

                end
                
                if fillPixels
                    % Perform Interpolation of non-viable pixels
                    % (nothing for now, design something to average the
                    % time series if possible
                    
                end

                mapDataClusters = NaN(flip(mapSize));
                idglobal = (mapIDglobal == j_dir);
                idxKtemp = idxK(idglobal);
                pixelLogtemp = pixelLogGlobal(idglobal,:);
                for k_cluster = 1:numel(idxKtemp)
                    mapDataClusters(pixelLogtemp(k_cluster,1),pixelLogtemp(k_cluster,2)) = idxKtemp(k_cluster);
                end
                
                if exist('XA','var') && exist('YA','var')
                    XPlot = XA./1e-6;
                    xlab = sprintf('X Position [\\mum]');
                    xlims = [0 max(scanSize)./1e-6];
                    YPlot = YA./1e-6;
                    ylab = sprintf('Y Position [\\mum]');
                    ylims = [0 max(scanSize)./1e-6];
                else
                    XPlot = X;
                    xlab = 'X Index';
                    xlims = [1 max(mapSize)];
                    YPlot = Y;
                    ylab = 'Y Index';
                    ylims = [1 max(mapSize)];
                end
                
                tiledlayout(n_rows,n_cols, 'padding', 'none', ...
                    'TileSpacing', 'compact', ...
                    'OuterPosition', [0 0.15 1 0.85])

                ax = nexttile;

                surf(XPlot,YPlot,rot90(heightImg,1),rot90(heightImg,1),'EdgeColor',mapEdgeCol)
                colormap(ax,'turbo')
                hold on
                title('Topography')
                ylabel(ylab)
                xlabel(xlab)
                xlim(xlims)
                ylim(ylims)
                cb = colorbar;
                caxis([0 climHeight]); % Absolute scale
                temp = (cb.Ticks' ./ 1e-6);
                for ii = 1:numel(temp)
                   cb.TickLabels{ii} = sprintf('%g \\mum',temp(ii));
                end
                view(2)
                pbaspect([1 1 1])
                hold off
                
                if plotIndentation

                    ax = nexttile;

                    surf(XPlot,YPlot,mapDataHeight,mapDataInd,'EdgeColor',mapEdgeCol)
                    colormap(ax,'turbo')
                    hold on
                    title(['Indentation' sprintf(', %d Hz',evalPt)])
                    ylabel(ylab)
                    xlabel(xlab)
                    xlim(xlims)
                    ylim(ylims)
                    cb = colorbar;
                    caxis([0 climInd]); % Absolute scale
                    temp = (cb.Ticks' ./ 1e-9);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%g nm',temp(ii));
                    end
                    view(2)
                    pbaspect([1 1 1])
                    hold off

                end

                ax = nexttile;

                switch clusterTarget
                    case 'force'
                        mapData = mapDataForce;
                        plotTitle = 'Observed Force Clustering';
                    case 'indentation'
                        mapData = mapDataInd;
                        plotTitle = 'Observed Indentation Clustering';
                    case 'storage'
                        mapData = mapDataStorage;
                        plotTitle = 'Storage Modulus Clustering';
                    case 'loss'
                        mapData = mapDataLoss;
                        plotTitle = 'Loss Modulus Clustering';
                    case 'angle'
                        mapData = mapDataAngle;
                        plotTitle = 'Loss Angle Clustering';
                    case 'relaxance'
                        mapData = mapDataRelaxance;
                        plotTitle = 'Relaxance Clustering';
                end
                
                surf(XPlot,YPlot,mapDataHeight,mapData,'EdgeColor',mapEdgeCol)
                colormap(ax,mapColorName)
                hold on
                title([plotTitle sprintf(', %d Hz',evalPt)])
                ylabel(ylab)
                xlabel(xlab)
                xlim(xlims)
                ylim(ylims)
                cb = colorbar;
                switch clusterTarget
                    case 'force'
                        caxis([0 10^ceil(log10(max(mapData,[],'all','omitnan')))]);
                        temp = (cb.Ticks' .* 1e-9);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%1.1g nN',temp(ii));
                        end
                    case 'indentation'
                        caxis([0 10^ceil(log10(max(mapData,[],'all','omitnan')))]);
                        temp = (cb.Ticks' .* 1e-9);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%1.1g nm',temp(ii));
                        end
                    case 'storage'
                        caxis([0 10^ceil(log10(max(mapData,[],'all','omitnan')))]);
                        temp = (cb.Ticks' .* 1e-3);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                        end
                    case 'loss'
                        caxis([0 10^ceil(log10(max(mapData,[],'all','omitnan')))]);
                        temp = (cb.Ticks' .* 1e-3);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                        end
                    case 'angle'
                        cb = colorbar;
                        cb.Ruler.TickLabelFormat='%d Deg';
                        caxis([0 90]);
                end
                                
                view(2)
                pbaspect([1 1 1])
                hold off
                
                ax = nexttile;

                surf(XPlot,YPlot,mapDataHeight,mapDataClusters,'EdgeColor',mapEdgeCol)
                colormap(ax,mapColorName)
                hold on
                title('DTW Clusters')
                ylabel(ylab)
                xlabel(xlab)
                xlim(xlims)
                ylim(ylims)
                if isempty(manualK)
                    cb = colorbar('Ticks',1:eva.OptimalK,...
                        'TickLabels',sprintfc('Bin %d',[1:eva.OptimalK]));
                else
                    cb = colorbar('Ticks',1:manualK,...
                        'TickLabels',sprintfc('Bin %d',[1:manualK]));
                end
                view(2)
                pbaspect([1 1 1])
                hold off
                
                saveas(mapPlotWindow,[plotFile '.fig'])
%                 saveas(mapPlotKMeans,[plotFile '.jpg'])
                print(mapPlotWindow,[plotFile '.png'],'-dpng','-r300');
                
                % Save Clusters to the Results File
                if ~isfield(resultsStruct.(varNames{j}), 'clusterData')
                    % Create a placeholder struct where we store all of the
                    % cluster results for each type of analysis. This is
                    % critical, since it means we can run studies on
                    % different clustering methods and compare results
                    % later.
                    temp = struct;
                    for jj = 1:6
                        switch jj
                            case 1
                                temp(jj).clusterVar = 'force';
                                
                            case 2
                                temp(jj).clusterVar = 'indentation';
                                
                            case 3
                                temp(jj).clusterVar = 'storage';
                                
                            case 4
                                temp(jj).clusterVar = 'loss';
                                
                            case 5
                                temp(jj).clusterVar = 'angle';
                                
                            case 6
                                temp(jj).clusterVar = 'relaxance';                                
                        end
                        temp(jj).clusterMap = {};
                        temp(jj).globalClusterMap = {};
                        temp(jj).clusterMap2D = {};
                        temp(jj).globalClusterMap2D = [];
                        temp(jj).lastUpdate = '';
                    end
                    
                    resultsStruct.(varNames{j}).clusterData = temp;
                    
                end
                
                cid = find(strcmp({resultsStruct.(varNames{j}).clusterData.clusterVar}, clusterTarget));
                resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap = num2cell(idxKtemp');
                resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap2D = mapDataClusters;
                resultsStruct.(varNames{j}).clusterData(cid).globalCentroidLocs = centroidLocs;
                resultsStruct.(varNames{j}).clusterData(cid).lastUpdate = datestr(now);
                
                if isfield(resultsStruct.(varNames{j}), 'trueBinsMap')
                    
                    % This is a test dataset where we have the "true" bins
                    % available. Quickly determine whether we have any
                    % incorrect values, and from that imply our accuracy.
                    
                    nBins = max(unique(idxKtemp),[],'all');
                    
                    for jj = 1:numel([resultsStruct.(varNames{j}).clusterData(:)])

                        % Loop through all of the observables used
                        if isempty(resultsStruct.(varNames{j}).clusterData(jj).globalClusterMap2D)
                            % Obviously we didn't analyze this observable,
                            % because there are no results!
                            continue;
                        end
                        
                        binNums = perms(1:nBins);
                    
                        for kk = 1:size(binNums,1)

                            % The clustered bins might not be in the correct
                            % "orientation". This would mean that the clustering
                            % could have correctly separated the regions, BUT it
                            % assigned the innermost bin as #1 instead of #3 for
                            % example. So we will loop through the combos and
                            % ensure we have the highest accuracy calculation
                            % possible.

                            % Cycle through the bin orientations
                            outputMap = resultsStruct.(varNames{j}).clusterData(jj).globalClusterMap2D;
                            mapDataClusters = NaN(size(outputMap));
                            for kbin = 1:size(binNums,2)
                                mapDataClusters(outputMap == kbin) = binNums(kk,kbin);
                            end

                            binDelta = (mapDataClusters ~= resultsStruct.(varNames{j}).trueBinsMap);
                            temp = 1 - (sum(binDelta,'all') / numel(binDelta));

                            if kk == 1
                            
                                outStruct.(dirLabel).clusterData(jj).globalClusterMap2D = resultsStruct.(varNames{j}).clusterData(jj).globalClusterMap2D;
                                outStruct.(dirLabel).clusterData(jj).globalClusterAccuracy = resultsStruct.(varNames{j}).clusterData(jj).globalClusterAccuracy;
                                outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);

                            elseif temp > outStruct.(dirLabel).clusterData(jj).clusterAccuracy/100

                                clusterAcc = temp;
                                outStruct.(dirLabel).clusterData(jj).globalClusterMap2D = mapDataClusters;
                                outStruct.(dirLabel).clusterData(jj).globalClusterAccuracy = 100*clusterAcc;
                                outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);
    %                             fprintf('\nThe Clustering Accuracy (%s) was %.2f%%\n',resultsStruct.(varNames{j}).clusterData(jj).clusterVar,100*clusterAcc);

                            end

                        end
                        
                        fprintf('\nThe Global Clustering Accuracy was %3.2f%% (%s)\n',100*clusterAcc, resultsStruct.(varNames{j}).clusterData(jj).clusterVar);

                    end
                                        
                end
                
                % Add our clusters to the file
                save([Files(j_dir).folder filesep Files(j_dir).name],'-struct','resultsStruct','-v7.3');
                                
            end
            
            clearvars resultsStruct mapSize
                        
        end
        
        fprintf('Saved Output Plots/Files!\n');
    
    catch ERROR
        
        fprintf('ERROR Clustering Directory #%d of %d\n',i_dir,length(Folders));
        fprintf('The identifier was:\n%s',ERROR.identifier);
        fprintf('Message:%s\n',ERROR.message);
        fprintf('Line Number:%d\n',ERROR.stack(end).line);
        fprintf('Skipping to next directory...\n');
        
    end

end

fprintf('Global Clustering Complete!\n');

% Clear Previous Parpool
if ~isempty(gcp('nocreate'))
   % Get the current pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end

end