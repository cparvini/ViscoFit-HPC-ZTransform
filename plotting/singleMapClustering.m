function [] = singleMapClustering(originalPath,N_workers,clusterTarget,varargin)
%SINGLEMAPCLUSTERING Perform K-Medoids Clustering of QI Map
%   This function takes in a path argument, the target variable for
%   clustering with k-medoids, and a variable number of arguments used to
%   correct/shift the map representations. Note that the only non-boolean
%   option for varargin is evalPt, which is the frequency (in Hz) to use
%   for plotting the frequency-dependent properties from the QI map. The
%   full spectrum of frequencies available for each pixel are used in the
%   clustering process, and evalPt is exclusively for visualization
%   purposes. This function has a pair, globalMapClustering(), which will
%   perform nearly the same analysis except it will combine many maps into
%   one large clustering dataset.
    
% User-Defined Settings
correctTilt = true;
hideSubstrate = true;
zeroSubstrate = true;
optimizeFlattening = false;
fillPixels = true;
logSteps = true;
plotIndentation = true;
evalPt = 1000;
n_reps = 10; % number of clustering replicates
maxK = 10; % Max number of cluster bins
if nargin > 1
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        correctTilt = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        hideSubstrate = varargin{i};
                    end
                case 3
                    if ~isempty(varargin{i})
                        zeroSubstrate = varargin{i};                        
                    end
                case 4
                    if ~isempty(varargin{i})
                        optimizeFlattening = varargin{i};                        
                    end
                case 5
                    if ~isempty(varargin{i})
                        fillPixels = varargin{i};                        
                    end
                case 6
                    if ~isempty(varargin{i})
                        logSteps = varargin{i};                        
                    end
                case 7
                    if ~isempty(varargin{i})
                        plotIndentation = varargin{i};                        
                    end
                case 8
                    if ~isempty(varargin{i})
                        evalPt = varargin{i};                        
                    end
                case 9
                    if ~isempty(varargin{i})
                        n_reps = varargin{i};                        
                    end
                case 10
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
climMax = 2e5; % Pa
climHeight = 15e-6; % meters, the JPK Nanowizard has a 15um piezo 
climInd = 1000e-9; % meters
trimHeight = 100e-9;
dFreq = 200; % Hz, step size between frames, if discrete
n_freqs = 10; % frames, number of frames per order of magnitude
n_datapoints = 10;

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

% Issue wrapper which allows script to continue if there is an empty
% directory, or other issue during processing.
try
    
    % Begin looping through the directories or files
    for i_dir = 1:length(Folders)
        path = Folders{i_dir};
        Files = dir([path filesep '*Results*zTransform*.mat']);

        if isempty(Files)
            error('The directory you selected does not contain a Z-Transform QI map. Please verify your FitResults file is in that directory and the filename contains "zTransform".');
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
                
                plotFile = [path filesep fileLabels{j_dir} saveLabel 'Clustering-' varNames{j}...
                    mapType];

                pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
                if correctTilt
                    temp = fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate,[],optimizeFlattening);
                    pixelHeightArray = cell2mat(temp);
                else
                    pixelHeightArray = cell2mat(pixelHeight_cell);
                end

                [minHeight,~] = min(pixelHeightArray);
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
                clusteringData = NaN(numel(mapDataHeight),numel(magList));
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

                        % Resample to known array of frequencies                        
                        mapDataStorage(idx_pixel) = interp1(freq,modelStorage,evalPt,'makima',...
                            NaN);
                        mapDataLoss(idx_pixel) = interp1(freq,modelLoss,evalPt,'makima',...
                            NaN);
                        mapDataAngle(idx_pixel) = interp1(freq,modelAngle,evalPt,'makima',...
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

                    switch clusterTarget
                        case 'force'
                            clusterInterp = F_t;
                        case 'indentation'
                            clusterInterp = h_t;
                        case 'storage'
                            clusterInterp = modelStorage;
                        case 'loss'
                            clusterInterp = modelLoss;
                        case 'angle'
                            clusterInterp = modelAngle;
                        case 'relaxance'
                            clusterInterp = modelRelaxance;
                    end

                    try
                        % Resample to known array of frequencies
                        obsOut = interp1(freq,clusterInterp,magList,'makima',...
                            NaN);
                        clusteringData(k_pixels,:) = obsOut;
                    catch
                        % Do nothing
                    end
                    
                    mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                    xc = xc + 1;

                end
                
                if fillPixels
                    % Perform Interpolation of non-viable pixels
                    % (nothing for now, design something to average the
                    % time series if possible
                    
                end

                opts = statset('UseParallel',parallelSet,...
                    'MaxIter',1e2,...
                    'Display','off');
                tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
                    'Distance',@dtwf,...
                    'Replicates',n_reps);
                ids = all(isnan(clusteringData),2); % Find excluded pixels
                clusteringData(ids,:) = [];
                pixelLog(ids,:) = [];
                
                % Try all of the cluster configurations. We start with 3
                % clusters to account for any substrate pixels that were
                % not trimmed. This is not ideal, but unavoidable with
                % automatic large-scale analysis.
                eva = evalclusters(clusteringData,tempfunc,'CalinskiHarabasz',...
                    'klist',((3-hideSubstrate):maxK));
                
                fprintf('\nOptimal Number of Bins: %d\n\n',eva.OptimalK);
                
                idxK = tempfunc(clusteringData,eva.OptimalK);
                
                mapDataClusters = NaN(flip(mapSize));
                for k_cluster = 1:numel(idxK)
                    mapDataClusters(pixelLog(k_cluster,1),pixelLog(k_cluster,2)) = idxK(k_cluster);
                end
                
                tiledlayout(n_rows,n_cols, 'padding', 'none', ...
                    'TileSpacing', 'compact', ...
                    'OuterPosition', [0 0.15 1 0.85])

                ax = nexttile;

                surf(X,Y,rot90(heightImg,1),rot90(heightImg,1),'EdgeColor','interp')
                colormap(ax,'turbo')
                hold on
                title('Topography')
                ylabel('Y Index')
                xlabel('X Index')
                xlim([1 max(mapSize)])
                ylim([1 max(mapSize)])
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

                    surf(X,Y,mapDataHeight,mapDataInd,'EdgeColor','interp')
                    colormap(ax,'turbo')
                    hold on
                    title(['Indentation' sprintf(', %d Hz',evalPt)])
                    ylabel('Y Index')
                    xlabel('X Index')
                    xlim([1 max(mapSize)])
                    ylim([1 max(mapSize)])
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
                
                surf(X,Y,mapDataHeight,mapData,'EdgeColor','interp')
                colormap(ax,mapColorName)
                hold on
                title([plotTitle sprintf(', %d Hz',evalPt)])
                ylabel('Y Index')
                xlabel('X Index')
                xlim([1 max(mapSize)])
                ylim([1 max(mapSize)])
                cb = colorbar;
                switch clusterTarget
                    case 'force'
                        caxis([0 climMax]);
                        temp = (cb.Ticks' .* 1e-9);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%0.2g nN',temp(ii));
                        end
                    case 'indentation'
                        caxis([0 climMax]);
                        temp = (cb.Ticks' .* 1e-9);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%0.2g nm',temp(ii));
                        end
                    case 'storage'
                        caxis([0 climMax]);
                        temp = (cb.Ticks' .* 1e-3);
                        for ii = 1:numel(temp)
                           cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                        end
                    case 'loss'
                        caxis([0 climMax]);
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

                surf(X,Y,mapDataHeight,mapDataClusters,'EdgeColor','interp')
                colormap(ax,mapColorName)
                hold on
                title('DTW Clusters')
                ylabel('Y Index')
                xlabel('X Index')
                xlim([1 max(mapSize)])
                ylim([1 max(mapSize)])
                cb = colorbar('Ticks',1:eva.OptimalK,...
                    'TickLabels',sprintfc('Bin %d',[1:eva.OptimalK]));
                view(2)
                pbaspect([1 1 1])
                hold off
                
                saveas(mapPlotWindow,[plotFile '.fig'])
%                 saveas(mapPlotKMeans,[plotFile '.jpg'])
                print(mapPlotWindow,[plotFile '.png'],'-dpng','-r300');

            end

        end

    end
    
catch ERROR
        
    fprintf('ERROR Animating Directory #%d of %d\n',i_dir,length(Folders));
    fprintf('The identifier was:\n%s',ERROR.identifier);
    fprintf('Message:%s\n',ERROR.message);
    fprintf('Skipping to next directory...\n');

end

% Clear Previous Parpool
if ~isempty(gcp('nocreate'))
   % Get the current pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end

end