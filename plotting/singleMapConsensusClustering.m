function [] = singleMapConsensusClustering(originalPath,N_workers,clusterTarget,varargin)
%SINGLEMAPCONSENSUSCLUSTERING Perform K-Medoids Consensus Clustering of QI
%Map
%   This function takes in a path argument, the target variable for
%   clustering with k-medoids, and a variable number of arguments used to
%   correct/shift the map representations. Note that the only non-boolean
%   option for varargin is evalPt, which is the frequency (in Hz) to use
%   for plotting the frequency-dependent properties from the QI map. The
%   full spectrum of frequencies available for each pixel are used in the
%   clustering process, and evalPt is exclusively for visualization
%   purposes.

% User-Defined Settings
hideSubstrate = true;
fillPixels = true;
logSteps = true;
plotIndentation = true;
evalPt = 1000;
n_reps = 100; % number of clustering replicates
maxK = 10; % Max number of cluster bins
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
n_freqs = 50; % frames, number of frames per order of magnitude
n_datapoints = 10;
clusterTargetList = {'force','indentation','storage','loss','relaxance'};

% Check to see if there are subdirectories
dirContents = dir(originalPath);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root filesep sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
else
    Folders = {originalPath};
end

% Define error functions
sse_global = @(data,model) sum((data-model).^2,'all');
mse_global = @(data,model,n) sum((data-model).^2,'all')./(length(data)-n);

% Clear Previous Parpool
if ~isempty(gcp('nocreate'))
   % Get the current pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end

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
    
% Begin looping through the directories or files
for i_dir = 1:length(Folders)
    
    % Issue wrapper which allows script to continue if there is an empty
    % directory, or other issue during processing.
    try
        
        path = Folders{i_dir};
        Files = dir([path filesep '*Results*zTransform*.mat']);

        if isempty(Files)
            error('The directory you selected does not contain a Z-Transform QI map. Please verify your results file is in that directory and the filename contains "zTransform".');
        end

        fileLabels = cell(numel(Files),1);
        for j = 1:numel(Files)
            temp = strsplit(Files(j).name,{'-','_','.'});
            idx = find(contains(lower(temp),{'fitresults','mapresults'}),1);
            fileLabels{j} = strjoin(temp([1 idx-1]),'-');
        end

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

                minFreq = Inf;
                maxFreq = 0;
                for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = min(tempf,[],'omitnan');
                    if any([isnan(temp),isempty(temp)]) || (numel(resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels})<n_datapoints)
                        continue;
                    end
                    if temp < minFreq
                        minFreq = temp;
                    end
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = max(tempf,[],'omitnan');
                    if isnan(temp)|| isempty(temp)
                        continue;
                    end
                    if temp > maxFreq
                        maxFreq = temp;
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
                            dFreq = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_freqs);
                        end

                        if temp >= tempmax
                            tempmax = temp*10;
                            dFreq = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_freqs);
                        end

                        tempf = 10.^( ( log10(temp) ) );
                        magList = horzcat(magList,tempf);
                        temp = temp + dFreq;
                    end
                end
                magList = unique(magList);
                freqList = flip(magList);
                timeList = 1./magList;

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
                
                saveLabel = '-Consensus';
                
                plotFile = [path filesep fileLabels{j_dir} saveLabel 'Clustering-' varNames{j}...
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
                clusteringData = NaN(numel(mapDataHeight),numel(magList),numel(clusterTargetList));
                pixelLog = NaN(numel(mapDataHeight),2);

                for k_pixels = 1:numel(heightImg)

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
                    
                    clusterInterp = [];
                    
                    for i_target = 1:numel(clusterTargetList)
                        
                        clusterTarget = clusterTargetList{i_target};
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
                            clusteringData(k_pixels,:,i_target) = obsOut;
                        catch
                            % Do nothing
                        end
                    
                    end
                        
                    mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                    xc = xc + 1;

                end
                
                if fillPixels
                    % Perform Interpolation of non-viable pixels
                    % (nothing for now, design something to average the
                    % time series if possible
                    
                end
                
                for k_loop = (3-hideSubstrate):maxK
                
                    % Pre-allocate consensus matrix
                    DCons = cell(numel(clusterTargetList),1);
                    idxKall = cell(size(DCons));
                    
                    for i_target = 1:numel(clusterTargetList)

                        clusteringDataLoop = clusteringData(:,:,i_target);
                        clusterTarget = clusterTargetList{i_target};

                        % Perform some data cleaning. Columns where we don't have enough
                        % observations compared to bins are a PROBLEM. We have to remove
                        % these before analyzing.
                        dataVerify = ~isnan(clusteringDataLoop);
                        goodCols = sum(dataVerify,1);
                        idRem = goodCols < maxK;

                        magListTemp = magList;
                        freqListTemp = freqList;
                        timeListTemp = timeList;
                        
                        magListTemp(idRem) = [];
                        freqListTemp(idRem) = [];
                        timeListTemp(idRem) = [];
                        clusteringDataLoop(:,idRem) = [];

                        opts = statset('UseParallel',parallelSet,...
                            'MaxIter',1e2,...
                            'Display','off');
                        tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
                            'Distance',@dtwf,...
                            'Replicates',1);
                        ids = all(isnan(clusteringDataLoop),2); % Find excluded pixels
                        clusteringDataLoop(ids,:) = [];
                        pixelLog(ids,:) = [];

                        idxK = NaN(size(clusteringDataLoop,1),n_reps);
                        
                        for i_rep = 1:n_reps
                        
                            % Pre-allocate
                            tempidx = NaN(size(clusteringDataLoop,1),1);
                            
                            % Perform clustering replicate
                            tempidx = tempfunc(clusteringDataLoop,k_loop);
                            
                            if i_rep > 1
                                
                                binNums = perms(1:k_loop);
                                clusterAcc = 0;
                                
                                for j_rep = 1:size(binNums,1)
                                    
                                    % Cycle through the bin orientations
                                    oldMap = idxK(:,i_rep-1);
                                    mapClusters = NaN(size(oldMap));
                                    for kbin = 1:size(binNums,2)
                                        mapClusters(oldMap == kbin) = binNums(j_rep,kbin);
                                    end

                                    binDelta = (mapClusters ~= oldMap);
                                    temp = 1 - (sum(binDelta,'all') / numel(binDelta));

                                    if j_rep == 1
                                        
                                        % Do nothing
                                        
                                    elseif temp > clusterAcc
                                        
                                        % The cluster assignments match the
                                        % previous map better. Overwrite
                                        % our old configuration
                                        clusterAcc = temp;
                                        tempidx = mapClusters;
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            idxK(:,i_rep) = tempidx;
                            
                        end
                        
                        % Create our consensus matrix for this k value
                        D = [];
                        D = squeeze(sum(bsxfun(@eq, idxK', permute(idxK', [1 3 2]))));
                        D = D/i_rep;
                        D(logical(eye(size(D)))) = 1;
                        
                        DCons{i_target} = D;
                        idxKall{i_target} = idxK;
                        
                    end
                    
                end
                
                % Find our optimal K using the elbow method
                
                
                % Create dissimilarity matrix based on the consensus matrix
                % for that run
                
                idxKout = [];
                
                % Visualize our results
                mapDataClusters = NaN(flip(mapSize));
                for k_cluster = 1:numel(idxK)
                    mapDataClusters(pixelLog(k_cluster,1),pixelLog(k_cluster,2)) = idxK(k_cluster);
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
                           cb.TickLabels{ii} = sprintf('%1.2g nN',temp(ii));
                        end
                    case 'indentation'
                        caxis([0 10^ceil(log10(max(mapData,[],'all','omitnan')))]);
                        temp = (cb.Ticks' .* 1e-9);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%1.2g nm',temp(ii));
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
                cb = colorbar('Ticks',1:eva.OptimalK,...
                    'TickLabels',sprintfc('Bin %d',[1:eva.OptimalK]));
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
                        temp(jj).clusterMap2D = [];
                        temp(jj).lastUpdate = '';
                    end
                    
                    resultsStruct.(varNames{j}).clusterData = temp;
                    
                end
                
                cid = find(strcmp({resultsStruct.(varNames{j}).clusterData.clusterVar}, clusterTarget));
                resultsStruct.(varNames{j}).clusterData(cid).clusterMap = num2cell(idxK');
                resultsStruct.(varNames{j}).clusterData(cid).clusterMap2D = mapDataClusters;
                resultsStruct.(varNames{j}).clusterData(cid).lastUpdate = datestr(now);
                
                if isfield(resultsStruct.(varNames{j}), 'trueBinsMap')
                    
                    % This is a test dataset where we have the "true" bins
                    % available. Quickly determine whether we have any
                    % incorrect values, and from that imply our accuracy.
                    
                    nBins = max(unique(idxK),[],'all');
                    
                    for jj = 1:numel([resultsStruct.(varNames{j}).clusterData(:)])

                        % Loop through all of the observables used
                        if isempty(resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D)
                            % Obviously we didn't analyze this observable,
                            % because there are no results!
                            continue;
                        end
                        
                        clusterAcc = 0;
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
                            outputMap = resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D;
                            mapDataClusters = NaN(size(outputMap));
                            for kbin = 1:size(binNums,2)
                                mapDataClusters(outputMap == kbin) = binNums(kk,kbin);
                            end

                            binDelta = (mapDataClusters ~= resultsStruct.(varNames{j}).trueBinsMap);
                            temp = 1 - (sum(binDelta,'all') / numel(binDelta));

                            if kk == 1

                                outStruct.(dirLabel).clusterData(jj).clusterMap2D = resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D;
                                outStruct.(dirLabel).clusterData(jj).clusterAccuracy = resultsStruct.(varNames{j}).clusterData(jj).clusterAccuracy;
                                outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);

                            elseif temp > outStruct.(dirLabel).clusterData(jj).clusterAccuracy/100

                                clusterAcc = temp;
                                outStruct.(dirLabel).clusterData(jj).clusterMap2D = mapDataClusters;
                                outStruct.(dirLabel).clusterData(jj).clusterAccuracy = 100*clusterAcc;
                                outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);
    %                             fprintf('\nThe Clustering Accuracy (%s) was %.2f%%\n',resultsStruct.(varNames{j}).clusterData(jj).clusterVar,100*clusterAcc);

                            end

                        end
                        
                        fprintf('\nThe Clustering Accuracy was %3.2f%% (%s)\n',100*clusterAcc, resultsStruct.(varNames{j}).clusterData(jj).clusterVar);

                    end
                                        
                end
                
                % Overwrite the old results
                save([Files(j_dir).folder filesep Files(j_dir).name],'-struct','resultsStruct','-v7.3');
                
            end

        end
    
    catch ERROR
        
        fprintf('ERROR Clustering Directory #%d of %d\n',i_dir,length(Folders));
        fprintf('The identifier was:\n%s',ERROR.identifier);
        fprintf('Message:%s\n',ERROR.message);
        fprintf('Line Number:%d\n',ERROR.stack(end).line);
        fprintf('Skipping to next directory...\n');

    end
        
end

% Clear Previous Parpool
if ~isempty(gcp('nocreate'))
   % Get the current pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end

end