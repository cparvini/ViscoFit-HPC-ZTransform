function [] = makeZTomography(originalPath,nPlanes,varargin)
%makeZTomography Create a 3D Tomographic Volume using Z-Transform Method
%   This function takes in a path argument and will subsequently create a
%   3D Volume showing the viscoelastic properties from a Z-Transform
%   Results file. The results file is created using either the
%   "analyze_map_zTransform" function or "fit_map_zTransform".

% User-Defined Settings
fillPixels = true;
showLabels = true;
climMax = 2e5;
if nargin > 1
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        fillPixels = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        showLabels = varargin{i};                        
                    end
                case 3
                    if ~isempty(varargin{i})
                        if isa(varargin{i},'numeric')
                            climMax = varargin{i};
                        else
                            climMax = 2e5; % Pa
                            warning('You provided a non-numeric input for the stiffness plot limit to makeZTomography(). Using the default value instead (2e5 Pa)');
                        end
                    end
                otherwise
                    fprintf('Passed additional parameters to makeZTomography() which were not used.');
            end
        end
    end
end

% Permanent Settings
errortype = 'sse';
plotModelAlways = false;
figX = 0;
figY = 0;
maxCol = 5;
maxwid = get(0,'screensize');
maxheight = maxwid(4);
maxwid = maxwid(3);
mapColorName = 'turbo';
boxSetting = 'on';
hardClim = true;
zMax = 15e-6; % meters, the JPK Nanowizard has a 15um piezo 
stiffMax = 10*climMax; % Pa
trimHeight = 100e-9;
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
    
% Begin looping through the directories or files
for i_dir = 1:length(Folders)

    % Issue wrapper which allows script to continue if there is an empty
    % directory, or other issue during processing.
    try
    
        fprintf('\nBeginning Tomographic Mapping for Directory #%d of %d\n',i_dir,length(Folders));
        
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
                
                % Grab some relevant settings
                correctTilt = resultsStruct.(varNames{j}).correctTilt;
                hideSubstrate = resultsStruct.(varNames{j}).hideSubstrate;
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
                dFreqAll = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),1);
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
                    dFreqAll(k_pixels) = median(gradient(tempf));
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
                dFreq = round(mean(dFreqAll,'omitnan'));
                while temp < maxFreq
                    tempf = 10.^( ( log10(temp) ) );
                    magList = horzcat(magList,tempf);
                    temp = temp + dFreq;
                end
                magList = unique(magList);
                timeList = flip(1./(magList));
                freqList = flip(magList);
                
                if isfield(resultsStruct.(varNames{j}),'bestParams') && plotModelAlways
                    plotModel = true;
                else
                    plotModel = false;
                end

                figWid = min([maxwid maxheight]);
                figHeight = figWid;

                if plotModel
                    mapType = '-Model';
                else
                    mapType = '-Raw';
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

                % Figure out the best limits for each subfigure
                forceMax = 0;
                indMax = 0;
                storMax = 0;
                lossMax = 0;
                relaxMax = 0;
                for k_pixels = 1:numel(pixelHeightArray)
                    
                    % Skip the hidden pixels
                    if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                        continue;
                    end
                    
                    % Update Force
                    if mean(resultsStruct.(varNames{j}).forceMap{k_pixels},'all','omitnan') > forceMax
                        forceMax = mean(resultsStruct.(varNames{j}).forceMap{k_pixels},'all','omitnan');
                    end
                    
                    % Update Indentation
                    if mean(resultsStruct.(varNames{j}).indMap{k_pixels},'all','omitnan') > indMax
                        indMax = mean(resultsStruct.(varNames{j}).indMap{k_pixels},'all','omitnan');
                    end
                    
                    % Update Relaxance
                    if mean(resultsStruct.(varNames{j}).relaxanceMap{k_pixels},'all','omitnan') > relaxMax
                        relaxMax = mean(resultsStruct.(varNames{j}).relaxanceMap{k_pixels},'all','omitnan');
                    end
                    
                    % Update Storage Mod
                    temp1 = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                    if mean(temp1,'all','omitnan') > storMax
                        storMax = mean(temp1,'all','omitnan');
                    end
                    
                    % Update Loss Mod
                    temp2 = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                    if mean(temp2,'all','omitnan') > lossMax
                        lossMax = mean(temp2,'all','omitnan');
                    end
                    
                end
                forceMax = 10^(ceil(log10(forceMax)));
                indMax = 10^(ceil(log10(indMax)));
                storMax = 10^(ceil(log10(storMax)));
                lossMax = 10^(ceil(log10(lossMax)));
                relaxMax = 10^(ceil(log10(relaxMax)));
                
                % Make blank map data
                mapDataStorage = NaN([flip(mapSize) numel(timeList)]);
                mapDataLoss = NaN([flip(mapSize) numel(timeList)]);
                mapDataAngle = NaN([flip(mapSize) numel(timeList)]);
                mapDataRelaxance = NaN([flip(mapSize) numel(timeList)]);
                mapDataError = NaN([flip(mapSize) numel(timeList)]);
                mapDataTerms = NaN([flip(mapSize) numel(timeList)]);
                mapDataHeight = NaN([flip(mapSize) numel(timeList)]);
                mapDataInd = NaN([flip(mapSize) numel(timeList)]);
                
                % Prep the movie
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
                
                progressString = sprintf('Populating Tomographic Dataset\nParallel Analysis Running...');
                hbar = parfor_progressbar(numel(timeList),progressString);
                warning('off');
                
                parfor k_t = 1:numel(timeList)
                    % Create temporary holders
                    mapDataStorageTemp = NaN(flip(mapSize));
                    mapDataLossTemp = NaN(flip(mapSize));
                    mapDataAngleTemp = NaN(flip(mapSize));
                    mapDataRelaxanceTemp = NaN(flip(mapSize));
                    mapDataErrorTemp = NaN(flip(mapSize));
                    mapDataTermsTemp = NaN(flip(mapSize));
                    mapDataHeightTemp = NaN(flip(mapSize));
                    mapDataIndTemp = NaN(flip(mapSize));
                    
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

                            if hideSubstrate && any(ismember(i_z,pixelSkip))
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

                            if any(cellfun(@isempty,dataIn(1:6))) || any(isnan(resultsStruct.(varNames{j}).frequencyMap{i_z}))
                                F_hz_all{i_z} = NaN;
                                h_hz_all{i_z} = NaN;
                                continue;
                            end

                            [~,~,F_hz_all{i_z},~,h_hz_all{i_z},~,~] = zTransformCurve(dataIn,'none',0.05,resultsStruct.(varNames{j}).ViscoClass.thinSample);

                        end

                    end
                    
                    for k_pixels = 1:numel(mapDataStorageTemp)

                        % Get the current pixel position
                        if xc > mapSize(1)
                            xc = 1;
                            yc = yc + 1;
                        end
                        
%                         idx_pixel = sub2ind([flip(mapSize) size(mapDataStorage,3)],mapSize(2)-yc,xc,k_t);
                        idx_pixel_slice = sub2ind(flip(mapSize),mapSize(2)-yc,xc);

                        if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                            xc = xc + 1;
                            continue;
                        end

                        if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                            xc = xc + 1;
                            continue;
                        end

                        if ~isfield(resultsStruct.(varNames{j}),'indMap')
                            F_hz = abs(F_hz_all{k_pixels});
                            h_hz = abs(h_hz_all{k_pixels});
                        else
                            F_hz = abs(resultsStruct.(varNames{j}).forceMap{k_pixels});
                            h_hz = abs(resultsStruct.(varNames{j}).indMap{k_pixels});
                        end

                        % Load and perform peak correction
                        freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                        [~,maxid] = max(F_hz);
                        freqAdj = freq(maxid);
                        freq = freq - freqAdj;

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

                            if (freqList(k_t) > max(freq,[],'omitnan')) || (freqList(k_t) < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq) < 2
                                mapDataStorageTemp(idx_pixel_slice) = NaN;
                                mapDataLossTemp(idx_pixel_slice) = NaN;
                                mapDataAngleTemp(idx_pixel_slice) = NaN;
                                mapDataRelaxanceTemp(idx_pixel_slice) = NaN;
                                mapDataErrorTemp(idx_pixel_slice) = NaN;
                                mapDataTermsTemp(idx_pixel_slice) = NaN;
                                mapDataHeightTemp(idx_pixel_slice) = NaN;
                                mapDataIndTemp(idx_pixel_slice) = NaN;
                                xc = xc + 1;
                            else
                                % Resample to known array of frequencies
                                if hideSubstrate && (max([interp1(freq,modelStorage,1./timeList(k_t),'makima',...
                                    NaN) interp1(freq,modelLoss,1./timeList(k_t),'makima',...
                                    NaN) 0],[],'omitnan') > stiffMax)

                                    xc = xc + 1;
                                    continue;
                                end

                                % Create our arrays in the time domain
                                t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                                h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                                mapDataIndTemp(idx_pixel_slice) = interp1(t_t,h_t,timeList(k_t),'makima',...
                                    1e-12);
                                
                                evalPt = 1/timeList(k_t);
                                mapDataStorageTemp(idx_pixel_slice) = interp1(freq,modelStorage,evalPt,'makima',...
                                    NaN);
                                mapDataLossTemp(idx_pixel_slice) = interp1(freq,modelLoss,evalPt,'makima',...
                                    NaN);
                                mapDataAngleTemp(idx_pixel_slice) = interp1(freq,modelAngle,evalPt,'makima',...
                                    NaN);
                                mapDataRelaxanceTemp(idx_pixel_slice) = interp1(freq,modelRelaxance,evalPt,'makima',...
                                    NaN);
                                
                                mapDataErrorTemp(idx_pixel_slice) = modelErrorTime;
                                mapDataTermsTemp(idx_pixel_slice) = bestidx;
                                mapDataHeightTemp(idx_pixel_slice) = pixelHeightArray(idx_pixel_slice);
                                xc = xc + 1;
                            end

                        else

                            % We are just plotting the data captured
                            % DIRECTLY from the z-transform method.
                            modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelAngle = atand(modelLoss./modelStorage);
                            modelRelaxance = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};

                            if (1/timeList(k_t) > max(freq,[],'omitnan')) || (1/timeList(k_t) < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq((freq >= min(freqList)) & (freq <= max(freqList)))) < 2
                                
                                mapDataStorageTemp(idx_pixel_slice) = NaN;
                                mapDataLossTemp(idx_pixel_slice) = NaN;
                                mapDataAngleTemp(idx_pixel_slice) = NaN;
                                mapDataRelaxanceTemp(idx_pixel_slice) = NaN;
                                mapDataHeightTemp(idx_pixel_slice) = NaN;
                                mapDataIndTemp(idx_pixel_slice) = NaN;
                                xc = xc + 1;
                                
                            else

                                % Resample to known array of frequencies
                                ids = ((freq >= min(freqList)) & (freq <= max(freqList)));

                                if hideSubstrate && (max([interp1(freq(ids),modelStorage(ids),1/timeList(k_t),'makima',...
                                    NaN) interp1(freq(ids),modelLoss(ids),1/timeList(k_t),'makima',...
                                    NaN) 0],[],'omitnan') > stiffMax)

                                    xc = xc + 1;
                                    continue;
                                end
                                
                                % Has to be from the time domain
                                % (normalization issues)
                                t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                                h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                                mapDataIndTemp(idx_pixel_slice) = interp1(t_t,h_t,timeList(k_t),'makima',...
                                    1e-12);
                                
                                evalPt = 1./timeList(k_t);
                                mapDataStorageTemp(idx_pixel_slice) = interp1(freq(ids),modelStorage(ids),evalPt,'makima',...
                                    NaN);
                                mapDataLossTemp(idx_pixel_slice) = interp1(freq(ids),modelLoss(ids),evalPt,'makima',...
                                    NaN);
                                mapDataAngleTemp(idx_pixel_slice) = interp1(freq(ids),modelAngle(ids),evalPt,'makima',...
                                    NaN);
                                mapDataRelaxanceTemp(idx_pixel_slice) = interp1(freq(ids),modelRelaxance(ids),evalPt,'makima',...
                                    NaN);

                                mapDataHeightTemp(idx_pixel_slice) = pixelHeightArray(idx_pixel_slice);
                                xc = xc + 1;
                            end

                        end

                        if fillPixels
                            % Perform Interpolation of non-viable pixels
                            % Start by creating masks for the pixels to fill
                            [m, n] = size(mapDataHeight(:,:,k_t));
                            neighbors4 = [-1, 1, m, -m];
                            neighbors8 = [neighbors4, -m-1, -m+1, m-1, m+1];
                            f4 = @(padimg,ind) mean(padimg(ind + neighbors4),"omitnan");
                            f8 = @(padimg,ind) mean(padimg(ind + neighbors8),"omitnan");

                            storageTemp = mapDataStorage(:,:,k_t);
                            storageTemp(pixelSkip) = 0;
                            paddedStorage = padarray(storageTemp, [1 1], 0, 'both');

                            lossTemp = mapDataLoss(:,:,k_t);
                            lossTemp(pixelSkip) = 0;
                            paddedLoss = padarray(lossTemp, [1 1], 0, 'both');

                            angleTemp = mapDataAngle(:,:,k_t);
                            angleTemp(pixelSkip) = 0;
                            paddedAngle = padarray(angleTemp, [1 1], 0, 'both');

                            % Storage Modulus Interpolation
                            originalNaNPos = find(isnan(storageTemp));
                            paddedImageNaNs = find(isnan(paddedStorage));
                            imglocav = storageTemp;

                            while nnz(isnan(imglocav)) > 0
                                originalNaNPos = find(isnan(imglocav));
                                for ii = 1:numel(originalNaNPos)
                                    imglocav(originalNaNPos(ii)) = f8(paddedStorage,paddedImageNaNs(ii));
                                end
                            end

                            storageTemp(originalNaNPos) = imglocav(originalNaNPos);
                            mapDataStorage(:,:,k_t) = storageTemp;

                            % Loss Modulus Interpolation
                            originalNaNPos = find(isnan(lossTemp));
                            paddedImageNaNs = find(isnan(paddedLoss));
                            imglocav = lossTemp;

                            while nnz(isnan(imglocav)) > 0
                                originalNaNPos = find(isnan(imglocav));
                                for ii = 1:numel(originalNaNPos)
                                    imglocav(originalNaNPos(ii)) = f8(paddedLoss,paddedImageNaNs(ii));
                                end
                            end

                            lossTemp(originalNaNPos) = imglocav(originalNaNPos);
                            mapDataLoss(:,:,k_t) = lossTemp;
                            
                            % Loss Angle Interpolation
                            originalNaNPos = find(isnan(angleTemp));
                            paddedImageNaNs = find(isnan(paddedAngle));
                            imglocav = angleTemp;

                            while nnz(isnan(imglocav)) > 0
                                originalNaNPos = find(isnan(imglocav));
                                for ii = 1:numel(originalNaNPos)
                                    imglocav(originalNaNPos(ii)) = f8(paddedAngle,paddedImageNaNs(ii));
                                end
                            end

                            angleTemp(originalNaNPos) = imglocav(originalNaNPos);
                            mapDataAngle(:,:,k_t) = angleTemp;
                            
                        end
                        
                    end
                    
                    % Data Output
                    mapDataStorage(:,:,k_t) = mapDataStorageTemp;
                    mapDataLoss(:,:,k_t) = mapDataLossTemp;
                    mapDataAngle(:,:,k_t) = mapDataAngleTemp;
                    mapDataRelaxance(:,:,k_t) = mapDataRelaxanceTemp;
                    mapDataError(:,:,k_t) = mapDataErrorTemp;
                    mapDataTerms(:,:,k_t) = mapDataTermsTemp;
                    mapDataHeight(:,:,k_t) = mapDataHeightTemp;
                    mapDataInd(:,:,k_t) = mapDataIndTemp;
                    
                    hbar.iterate(1) % Increase progressbar by 1 iteration
                end
                
                close(hbar)
                warning('on');
                
                mapDataStorage(mapDataStorage == 0) = NaN;
                mapDataLoss(mapDataLoss == 0) = NaN;
                mapDataAngle(mapDataAngle == 0) = NaN;
                mapDataRelaxance(mapDataRelaxance == 0) = NaN;
                mapDataInd(mapDataInd == 0) = NaN;
                
                if exist('XA','var') && exist('YA','var')
                    XPlot = XA./1e-6;
                    xlab = sprintf('X Position [\\mum]');
                    xlims = [0 max(scanSize)./1e-6];
                    YPlot = YA./1e-6;
                    ylab = sprintf('Y Position [\\mum]');
                    ylims = [0 max(scanSize)./1e-6];
                    zlab = sprintf('Z Position [\\mum]');
                    zlims = [0 zMax./1e-6];
                    zscale = 1e6;
                else
                    XPlot = X;
                    xlab = 'X Index';
                    xlims = [1 max(mapSize)];
                    YPlot = Y;
                    ylab = 'Y Index';
                    ylims = [1 max(mapSize)];
                    zlab = 'Z Index';
                    zlims = [0 zMax./1e-6];
                    zscale = 1;
                end
                
                % Find our slice locations
                switch nPlanes
                    case 0
                        xslice = [];
                        yslice = [];
                        zslice = [];
                        
                    case 1
                        xslice = (xlims(2)-xlims(1))/2;
                        yslice = (ylims(2)-ylims(1))/2;
                        zslice = [];
                        
                    otherwise
                        xslice = linspace(xlims(1),xlims(2),nPlanes);
                        yslice = linspace(ylims(1),ylims(2),nPlanes);
                        zslice = [];
                end

                for i_plot = 1:4

                    figure(mapPlotWindow)
                    clf
                    
                    switch i_plot
                        case 1
                            mapData = mapDataStorage;
                            plotTitle = 'Storage Modulus';
                            maxScale = storMax;
                            saveLabel = 'Storage';
                            
                        case 2
                            mapData = mapDataLoss;
                            plotTitle = 'Loss Modulus';
                            maxScale = lossMax;
                            saveLabel = 'Loss';
                            
                        case 3
                            mapData = mapDataAngle;
                            plotTitle = 'Loss Angle';
                            maxScale = 90;
                            saveLabel = 'Angle';
                            
                        case 4
                            mapData = abs(mapDataRelaxance);
                            plotTitle = 'Relaxance';
                            maxScale = relaxMax;
                            saveLabel = 'Relaxance';
                            
                    end

                    for j_plot = 1:size(mapData,3)
                        
                        tempX = XPlot;
                        tempY = YPlot;
                        tempZ = zscale*(mapDataHeight(:,:,j_plot)-mapDataInd(:,:,j_plot));
                        tempV = mapData(:,:,j_plot);
                        
                        xbest = NaN(nPlanes,1);
                        ybest = NaN(nPlanes,1);
                        for k_plot = 1:nPlanes
                            [~,xbest(k_plot)] = min(abs(tempX(1,:)-xslice(k_plot)));
                            [~,ybest(k_plot)] = min(abs(tempY(:,1)-yslice(k_plot)));
                        end
                        idslice = (ismember(tempX,tempX(1,xbest))) | (ismember(tempY,tempY(ybest,1)));
                        tempX(~idslice) = NaN;
                        tempY(~idslice) = NaN;
                        tempZ(~idslice) = NaN;
                        tempV(~idslice) = NaN;
                        
                        hold on
%                         surf(XPlot(:,:,j_plot),YPlot(:,:,j_plot),mapDataHeight(:,:,j_plot)-mapDataInd(:,:,j_plot),mapData(:,:,j_plot),'EdgeColor','interp')
                        surf(tempX,tempY,tempZ,tempV,'EdgeColor','interp','FaceColor','none')
                        hold off
                        
                    end
                    
                    hold on
                    surf(XPlot,YPlot,zscale*mapDataHeight(:,:,1),max(mapData,[],3),'FaceAlpha',0.3,'EdgeAlpha',0.3)
                    hold off
                    
                    ax = gca;
                    colormap(ax,mapColorName)
                    box(ax,boxSetting)
                    grid minor
                    title(plotTitle)
                    if showLabels
                        xlabel(xlab)
                        ylabel(ylab)
                        zlabel(zlab)
                    else
                        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
                    end
                    xlim(xlims)
                    ylim(ylims)
                    zlim(zlims)
                    cb = colorbar;
                    
                    switch i_plot
                        case 1
                            if ~hardClim
                                caxis([0 maxScale]); % Absolute scale
                            else
                                caxis([0 climMax]);
                            end
                            temp = (cb.Ticks' .* 1e-3);
                            for ii = 1:numel(temp)
                               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                            end
                            
                        case 2
                            if ~hardClim
                                caxis([0 maxScale]); % Absolute scale
                            else
                                caxis([0 climMax]);
                            end
                            temp = (cb.Ticks' .* 1e-3);
                            for ii = 1:numel(temp)
                               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                            end
                            
                        case 3
                            cb.Ruler.TickLabelFormat='%g Deg';
                            caxis([0 90]);
                            
                        case 4
                            if ~hardClim
                                caxis([0 maxScale]); % Absolute scale
                            else
                                caxis([0 climMax]);
                            end
                            temp = (cb.Ticks' .* 1e-3);
                            for ii = 1:numel(temp)
                               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                            end
                            
                    end
                    
                    view(3)
                    pbaspect([1 1 1])
                    drawnow
                    
                    % Create our filename
                    plotFile = [path filesep fileLabels{j_dir} '-MapTomography-' varNames{j}...
                        mapType '-' saveLabel];
                    
                    % Using Print for Higher Quality
                    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
                    exportgraphics(ax,[plotFile '.png'],'Resolution',300);

                end

                if plotModel
                    
                    figure(mapPlotWindow)
                    clf
                    
                    mapData = mapDataError;
                    plotTitle = 'Model Error';
                    saveLabel = 'ModelError';
                    
                    for j_plot = 1:size(mapData,3)
                        
                        tempX = XPlot;
                        tempY = YPlot;
                        tempZ = zscale*(mapDataHeight(:,:,j_plot)-mapDataInd(:,:,j_plot));
                        tempV = mapData(:,:,j_plot);
                        
                        xbest = NaN(nPlanes,1);
                        ybest = NaN(nPlanes,1);
                        for k_plot = 1:nPlanes
                            [~,xbest(k_plot)] = min(abs(tempX(1,:)-xslice(k_plot)));
                            [~,ybest(k_plot)] = min(abs(tempY(:,1)-yslice(k_plot)));
                        end
                        idslice = (ismember(tempX,tempX(1,xbest))) | (ismember(tempY,tempY(ybest,1)));
                        tempX(~idslice) = NaN;
                        tempY(~idslice) = NaN;
                        tempZ(~idslice) = NaN;
                        tempV(~idslice) = NaN;
                        
                        hold on
%                         surf(XPlot(:,:,j_plot),YPlot(:,:,j_plot),mapDataHeight(:,:,j_plot)-mapDataInd(:,:,j_plot),mapData(:,:,j_plot),'EdgeColor','interp')
                        surf(tempX,tempY,tempZ,tempV,'EdgeColor','interp','FaceColor','none')
                        hold off
                        
                    end
                    
                    hold on
                    surf(XPlot,YPlot,zscale*mapDataHeight(:,:,1),max(mapData,[],3),'FaceAlpha',0.3,'EdgeAlpha',0.3)
                    hold off
                    
                    ax = gca;
                    colormap(ax,mapColorName)
                    box(ax,boxSetting)
                    title(plotTitle)
                    if showLabels
                        xlabel(xlab)
                        ylabel(ylab)
                        zlabel(zlab)
                    else
                        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[])
                    end
                    xlim(xlims)
                    ylim(ylims)
                    zlim(zlims)
                    colorbar
                    
                    view(3)
                    pbaspect([1 1 1])
                    drawnow
                    
                    % Create our filename
                    plotFile = [path filesep fileLabels{j_dir} '-MapTomography-' varNames{j}...
                        mapType '-' saveLabel];
                    
                    % Using Print for Higher Quality
                    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
                    exportgraphics(ax,[plotFile '.png'],'Resolution',300);
                    
                end
                
            end

        end
        
        fprintf('Tomographic Mapping Complete for Directory #%d of %d\n',i_dir,length(Folders));

    catch ERROR
        
        fprintf('ERROR Creating Tomographic Maps for Directory #%d of %d\n',i_dir,length(Folders));
        fprintf('The identifier was:\n%s',ERROR.identifier);
        fprintf('Message:%s\n',ERROR.message);
        fprintf('Line Number:%d\n',ERROR.stack(end).line);
        fprintf('Skipping to next directory...\n');    
    
    end

end

end

