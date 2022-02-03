clear all
close all
clc

addpath(genpath(['..' filesep 'lib']));
addpath(genpath(['..' filesep 'plotting']));

stillRunning = true;
while stillRunning
    
    % User-Defined Settings
    errortype = 'sse';
    hideSubstrate = true;
    zeroSubstrate = true; % In addition to correcting for tilt, also set the new, flat surface to have a minimum value starting at zero.
    fillPixels = true;
    plotIndentation = true;
    logSteps = true;
    figX = 0;
    figY = 50;
    nTicks = 5;
    plotModel = false;
    maxCol = 5;
    n_plots = 4+plotModel+plotIndentation;
    n_rows = ceil(n_plots/maxCol);
    n_cols = min([n_plots maxCol]);
    maxwid = get(0,'screensize');
    maxwid = maxwid(3);
    mult = min([400 maxwid/n_cols]);
    figWid = mult*n_cols;
    figHeight = max([mult*n_rows figWid/n_plots]);
    mapColorName = 'turbo';
    climMax = 2e5; % Pa
    climHeight = 15e-6; % meters
    climInd = 500e-9; % meters
    stiffMax = 10*climMax; % Pa
    trimHeight = 50e-9;
    dFreq = 200; % Hz, step size between frames
    n_frames = 150; % frames, number of frames per order of magnitude
    n_datapoints = 10;
    fps = 15;
    
    if plotModel
        mapType = '-Model';
    else
        mapType = '-Raw';
    end
    
    if exist('originalPath','var')
        startdir = originalPath;
    else
        startdir = pwd;
    end
    
    % Pick the AFM Data Directory and Choose Data Extraction Settings
    originalPath = uigetdir(startdir,...
            'Select the Folder Containing Your AFM Files');

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

    % Clear old figures if they exist
    if ~exist('mapPlotWindow','var')
        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
    else
        try
            figure(mapPlotWindow)
            clf
        catch
            clearvars mapPlotWindow
            mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
        end
    end

    % Make a list of unique colors for the test conditions
    directoryColors = num2cell(jet(numel(Folders)),2);
    
    % Define error functions
    sse_global = @(data,model) sum((data-model).^2,'all');
    mse_global = @(data,model,n) sum((data-model).^2,'all')./(length(data)-n);
    
    % Labels to search for
    varSearches = {'maxwellFit','voigtFit','PLRFit'};

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

            markerStyles = {'o','d','s'};
            colorStyles = {'r','b','g','m'};
            
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
                        if temp >= tempmax
                            tempmax = temp*10;
                            dFreq = round((tempmax-(tempmax/10))/n_frames);
                        end
                        
                        tempf = 10.^( ( log10(temp) ) );
                        magList = horzcat(magList,tempf);
                        temp = temp + dFreq;
                    end
                end
                magList = unique(magList);
                
%                 % Cleverly make a list of numbers to sample at which
%                 % divide each order of magnitude into 10 samples.
%                 freqList = [];
%                 for k = 2:numel(magList)
%                     freqList = horzcat(freqList,linspace(magList(k-1),magList(k),nFrames));
%                 end
%                 freqList = flip(unique(freqList)); % Process from high to low freq
                freqList = flip(magList);

                % Prep the movie
                u = uicontrol(mapPlotWindow,'Style','slider');
                u.Position = [20 15 figWid-230 20];
                u.Max = max(freqList);
                u.Min = min(freqList);
                u.Value = freqList(1);
                u2 = uicontrol(mapPlotWindow,'Style','edit');
                u2.Position = [figWid-175 10 150 40];
                u2.String = [num2str(round(freqList(1))) ' Hz'];
                u2.FontSize = 16;
                drawnow
                
                gifFile = [path filesep fileLabels{j_dir} '-MapAnimation-' varNames{j}...
                    mapType '.gif'];
                movieFile = [path filesep fileLabels{j_dir} '-MapMovie-' varNames{j}...
                    mapType '.mp4'];
                
                v = VideoWriter(movieFile,'MPEG-4');
                v.FrameRate = fps;
                v.Quality = 100;
                open(v);
                
                pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
                pixelHeightArray = cell2mat(fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate));

                [minHeight,~] = min(pixelHeightArray);
                substrateCutoff = minHeight + trimHeight;
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel(pixelHeightArray);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list
                
                heightImg = zeros(mapSize);
                
                for k_freq = 1:numel(freqList)

                    % Make blank map data
                    mapDataStorage = NaN(mapSize);
                    mapDataLoss = NaN(mapSize);
                    mapDataAngle = NaN(mapSize);
                    mapDataError = NaN(mapSize);
                    mapDataTerms = NaN(mapSize);
                    mapDataHeight = NaN(mapSize);
                    mapDataInd = NaN(mapSize);

                    % Position for the map
                    xc = 1;
                    yc = 0;

                    for k_pixels = 1:numel(mapDataStorage)
                        
                        % Get the current pixel position
                        if xc > mapSize(1)
                            xc = 1;
                            yc = yc + 1;
                        end
                        idx_pixel = sub2ind(mapSize,mapSize(2)-yc,xc);
                        
                        heightImg(idx_pixel) = pixelHeightArray(idx_pixel);
                        
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
                            
%                             [~,idx] = min(abs(freq-freqList(k_freq)));
                            modelRelaxance = modelStorage + 1j*modelLoss;
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq) < 2
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataRelaxance(idx_pixel) = NaN;
                                mapDataError(idx_pixel) = NaN;
                                mapDataTerms(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                mapDataInd(idx_pixel) = NaN;
                                xc = xc + 1;
                            else
                                % Resample to known array of frequencies
%                                 ids = ((freq >= min(freqList)) & (freq <= max(freqList)));

                                if hideSubstrate && (max([interp1(freq,modelStorage,freqList(k_freq),'makima',...
                                    NaN) interp1(freq,modelLoss,freqList(k_freq),'makima',...
                                    NaN) 0],[],'omitnan') > stiffMax)
                                    
                                    xc = xc + 1;
                                    continue;
                                end
                                
                                mapDataStorage(idx_pixel) = interp1(freq,modelStorage,freqList(k_freq),'makima',...
                                    NaN);
                                mapDataLoss(idx_pixel) = interp1(freq,modelLoss,freqList(k_freq),'makima',...
                                    NaN);
                                mapDataAngle(idx_pixel) = interp1(freq,modelAngle,freqList(k_freq),'makima',...
                                    NaN);
                                mapDataRelaxance(idx_pixel) = interp1(freq,modelRelaxance,freqList(k_freq),'makima',...
                                    NaN);
                                
                                % Has to be from the time domain
                                % (normalization issues)
                                evalPt = 1./(freqList(k_freq));
                                t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                                h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                                mapDataInd(idx_pixel) = interp1(t_t,h_t,evalPt,'makima',...
                                    1e-12);
                                                                
%                                 mapDataStorage(idx_pixel) = modelStorage(idx);
%                                 mapDataLoss(idx_pixel) = modelLoss(idx);
%                                 mapDataAngle(idx_pixel) = modelAngle(idx);
%                                 mapDataRelaxance(idx_pixel) = modelRelaxance(idx);
                                
                                mapDataError(idx_pixel) = modelErrorTime;
                                mapDataTerms(idx_pixel) = bestidx;
                                mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                                xc = xc + 1;
                            end
                            
                        else
                            
                            % We are just plotting the data captured
                            % DIRECTLY from the z-transform method.
                            modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelAngle = atand(modelLoss./modelStorage);
                                                        
%                             if numel(freq) > 1
%                                 % Smooth data
%                                 modelStorage = smoothdata(modelStorage,'gaussian',20);
%                                 modelLoss = smoothdata(modelLoss,'gaussian',20);
%                                 modelAngle = atand(modelLoss./modelStorage);
% 
%                                 figure
%                                 plot(freq,smoothStorage,'r',freq,modelStorage,'bo')
%                                 plot(freq,smoothLoss,'r',freq,modelLoss,'bo')
%                                 plot(freq,smoothAngle,'r',freq,modelAngle,'bo')
%                             end
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq((freq >= min(freqList)) & (freq <= max(freqList)))) < 2
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                mapDataInd(idx_pixel) = NaN;
                                xc = xc + 1;
                            else
                                
%                                 [~,idx] = min(abs(freq-freqList(k_freq)));                                
                                % Resample to known array of frequencies
                                ids = ((freq >= min(freqList)) & (freq <= max(freqList)));
                                
                                if hideSubstrate && (max([interp1(freq(ids),modelStorage(ids),freqList(k_freq),'makima',...
                                    NaN) interp1(freq(ids),modelLoss(ids),freqList(k_freq),'makima',...
                                    NaN) 0],[],'omitnan') > stiffMax)
                                    
                                    xc = xc + 1;
                                    continue;
                                end
                                
                                mapDataStorage(idx_pixel) = interp1(freq(ids),modelStorage(ids),freqList(k_freq),'makima',...
                                    NaN);
                                mapDataLoss(idx_pixel) = interp1(freq(ids),modelLoss(ids),freqList(k_freq),'makima',...
                                    NaN);
                                mapDataAngle(idx_pixel) = interp1(freq(ids),modelAngle(ids),freqList(k_freq),'makima',...
                                    NaN);
                                
                                % Has to be from the time domain
                                % (normalization issues)
                                evalPt = 1./(freqList(k_freq));
                                t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                                h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                                mapDataInd(idx_pixel) = interp1(t_t,h_t,evalPt,'makima',...
                                    1e-12);
                                
%                                 mapDataStorage(idx_pixel) = modelStorage(idx);
%                                 mapDataLoss(idx_pixel) = modelLoss(idx);
%                                 mapDataAngle(idx_pixel) = modelAngle(idx);

                                mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                                xc = xc + 1;
                            end
                        
                        end

                    end

                    if fillPixels
                        % Perform Interpolation of non-viable pixels
                        % Start by creating masks for the pixels to fill
                        [m, n] = size(mapDataHeight);
                        neighbors4 = [-1, 1, m, -m];
                        neighbors8 = [neighbors4, -m-1, -m+1, m-1, m+1];
                        f4 = @(padimg,ind) mean(padimg(ind + neighbors4),"omitnan");
                        f8 = @(padimg,ind) mean(padimg(ind + neighbors8),"omitnan");

                        storageTemp = mapDataStorage;
                        storageTemp(pixelSkip) = 0;
                        paddedStorage = padarray(storageTemp, [1 1], 0, 'both');
                        
                        lossTemp = mapDataLoss;
                        lossTemp(pixelSkip) = 0;
                        paddedLoss = padarray(lossTemp, [1 1], 0, 'both');
                        
                        angleTemp = mapDataAngle;
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
                        
                        mapDataStorage(originalNaNPos) = imglocav(originalNaNPos);
                        
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
                        
                        mapDataLoss(originalNaNPos) = imglocav(originalNaNPos);
                        
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
                        
                        mapDataAngle(originalNaNPos) = imglocav(originalNaNPos);
                        
%                         % Apply second method for further improvement
%                         storageNaNMask = mapDataStorage;
%                         storageNaNMask(pixelSkip) = 0; % Don't interpolate substrate pixels, if they exist
%                         storageNaNMask = imbinarize(storageNaNMask,"adaptive");
%                         storageNaNMask = bwareaopen(storageNaNMask,1);
% 
%                         lossNaNMask = mapDataLoss;
%                         lossNaNMask(pixelSkip) = 0; % Don't interpolate substrate pixels, if they exist
%                         lossNaNMask = imbinarize(lossNaNMask,"adaptive");
%                         lossNaNMask = bwareaopen(lossNaNMask,1);
% 
%                         angleNaNMask = mapDataAngle;
%                         angleNaNMask(pixelSkip) = 0; % Don't interpolate substrate pixels, if they exist
%                         angleNaNMask = imbinarize(angleNaNMask,"adaptive");
%                         angleNaNMask = bwareaopen(angleNaNMask,1);
% 
%                         % First, do the storage map
%                         edgeMask = edge(storageNaNMask);
%                         bwl = bwlabel(edgeMask);
%                         bwl2 = bwlabel(storageNaNMask);
%                         imX = mapDataStorage;
%                         imX(storageNaNMask) = NaN;
% 
%                         for ii = 1:max(bwl(:))
%                             thisEdgeMean = mean(imX(bwl == ii), "omitnan");
%                             imX(bwl2 == ii) = thisEdgeMean;
%                         end
%                         
%                         mapDataStorage(~storageNaNMask) = imX(~storageNaNMask);
%                         
%                         % Next, do the loss map
%                         edgeMask = edge(lossNaNMask);
%                         bwl = bwlabel(edgeMask);
%                         bwl2 = bwlabel(lossNaNMask);
%                         imX = mapDataLoss;
%                         imX(lossNaNMask) = NaN;
% 
%                         for ii = 1:max(bwl(:))
%                             thisEdgeMean = mean(imX(bwl == ii), "omitnan");
%                             imX(bwl2 == ii) = thisEdgeMean;
%                         end
%                         
%                         mapDataLoss(~lossNaNMask) = imX(~lossNaNMask);
%                         
%                         % Last, do the Angle map
%                         edgeMask = edge(angleNaNMask);
%                         bwl = bwlabel(edgeMask);
%                         bwl2 = bwlabel(angleNaNMask);
%                         imX = mapDataAngle;
%                         imX(angleNaNMask) = NaN;
% 
%                         for ii = 1:max(bwl(:))
%                             thisEdgeMean = mean(imX(bwl == ii), "omitnan");
%                             imX(bwl2 == ii) = thisEdgeMean;
%                         end
%                         
%                         mapDataAngle(~angleNaNMask) = imX(~angleNaNMask);

                    end
                    
                    mapDataStorage(mapDataStorage == 0) = NaN;
                    mapDataLoss(mapDataLoss == 0) = NaN;
                    mapDataAngle(mapDataAngle == 0) = NaN;
                    mapDataInd(mapDataInd == 0) = NaN;
                    
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
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
                    cb = colorbar;
%                     caxis([0 max(rot90(heightImg,1),[],'all')]); % Relative scale
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
                        title('Indentation')
                        ylabel('Y Index')
                        xlabel('X Index')
                        xlim([1 mapSize(1)])
                        ylim([1 mapSize(2)])
                        cb = colorbar;
%                         caxis([0 max(mapDataInd,[],'all')]); % Relative scale
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
                    
                    surf(X,Y,mapDataHeight,mapDataStorage,'EdgeColor','interp')
                    colormap(ax,mapColorName)
                    hold on
                    title('Storage Modulus')
                    ylabel('Y Index')
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataStorage,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    caxis([0 climMax]);
                    temp = (cb.Ticks' .* 1e-3);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                    end
                    view(2)
                    pbaspect([1 1 1])
                    hold off

                    ax = nexttile;
                    
                    surf(X,Y,mapDataHeight,mapDataLoss,'EdgeColor','interp')
                    colormap(ax,mapColorName)
                    hold on
                    title('Loss Modulus')
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataLoss,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    caxis([0 climMax]);
                    temp = (cb.Ticks' .* 1e-3);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
                    end
                    view(2)
                    pbaspect([1 1 1])
                    hold off
                    
                    ax = nexttile;
                    
                    surf(X,Y,mapDataHeight,mapDataAngle,'EdgeColor','interp')
                    colormap(ax,mapColorName)
                    hold on
                    title('Loss Angle')
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataAngle,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    cb.Ruler.TickLabelFormat='%g Deg';
                    caxis([0 90]);
                    view(2)
                    pbaspect([1 1 1])
                    hold off

                    if plotModel
                        ax = nexttile;
                        
                        surf(X,Y,mapDataHeight,mapDataTerms)
                        colormap(ax,mapColorName)
                        colormap(gca,'parula')
                        hold on
                        xlim([1 mapSize(1)])
                        ylim([1 mapSize(2)])
                        title(sprintf('Number of Terms'))
                        xlabel('X Index')
                        cb = colorbar;
                        set(cb,'YTick',1:numel(resultsStruct.(varNames{j}).bestParams))
                        view(2)
                        pbaspect([1 1 1])
                        hold off
                    end
                    
                    % Save Animation
                    u.Value = freqList(k_freq);
                    u2.String = [num2str(round(freqList(k_freq))) ' Hz'];
                    drawnow
                    M(k_freq) = getframe(mapPlotWindow);
                    frame = M(k_freq);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);

                    if k_freq == 1
                        imwrite(imind,cm,gifFile,'gif','DelayTime',0.05,'Loopcount',inf);
                    else
                        imwrite(imind,cm,gifFile,'gif','DelayTime',0.05,'WriteMode','append');
                    end
                    
                    % Write to mp4
                    writeVideo(v,frame);
                    
                end
                
                % End mp4
                close(v);
                
%                 % Display Movie
%                 playMov = true;
%                 qst = "Would you like to play the movie again?";
%                 while playMov
%                     movie(M);
%                     resp = questdlg(qst,"Replay Request",'Yes','No','No');
%                     if strcmp(resp,'No')
%                         playMov = false;
%                     end
%                 end

            end
            
        end
        
    end
    
    % Prompt user
    answer = questdlg('Would you like to analyze another directory?', ...
        'Options', ...
        'No','Yes','Yes');

    % Handle response
    switch answer
        case 'No'
            clearvars -except originalPath
            stillRunning = false;
        case 'Yes'
            close all
            clearvars -except stillRunning originalPath
    end
    
end
    
% Open the originally requested directory
winopen(originalPath);