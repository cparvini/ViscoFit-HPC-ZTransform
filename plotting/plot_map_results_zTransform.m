clear all
close all
clc

addpath(genpath('..\lib'));
addpath(genpath('..\plotting'));

stillRunning = true;
while stillRunning
    
    % User-Defined Settings
    nCols = 2;
    errortype = 'sse';
    figX = 100;
    figY = 100;
    nTicks = 5;
    plotModel = false;
    plotKMeans = true;
    gradientOrder = 2; % Number of derivatives to take of the relaxance map 
    kmeansBins = 5;
    lowfreqCutoff = sqrt(eps);
    figWid = 600*(3+plotModel);
    figHeight = figWid/(3+plotModel);
    mapColorName = 'parula';
    climMax = 6e5; % Pa
    
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
        Folders = cellfun(@(root,sub)[root '\' sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
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
    if ~exist('mapPlotKMeans','var')
        mapPlotKMeans = figure('Position',[figX figY figHeight figHeight]);
    else
        try
            figure(mapPlotKMeans)
            clf
        catch
            clearvars mapPlotKMeans
            mapPlotKMeans = figure('Position',[figX+100 figY figHeight figHeight]);
        end
    end
    if ~exist('mapPlotGradients','var')
        mapPlotGradients = figure('Position',[figX figY figHeight figHeight]);
    else
        try
            figure(mapPlotGradients)
            clf
        catch
            clearvars mapPlotGradients
            mapPlotGradients = figure('Position',[figX+100 figY figHeight figHeight]);
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
        Files = dir([path '\*Results*zTransform*.mat']);
        
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
            resultsStruct = load([Files(j_dir).folder '\' Files(j_dir).name],'-mat');

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
                
%                 freqList = [];
%                 for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
%                     temp = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
%                     temp = 10.^(unique(floor(log10(temp))));
%                     if any(isnan(temp)) || isempty(temp)
%                         continue;
%                     end
%                     
%                     % Remove frequencies to query where we have a curve
%                     % that doesn't contain enough info (so all pixels can
%                     % have data for our plots).
%                     if k_pixels == 1
%                         freqList = temp;
%                     elseif numel(temp) < numel(freqList)
%                         freqList = temp;
%                     end
%                 end
%                 
%                 freqList = freqList(freqList>0);
                
                minFreq = 0;
                maxFreq = Inf;
                numFreqPoints = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),1);
                minFreqPoints = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),1);
                maxFreqPoints = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),1);
                for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>lowfreqCutoff);
                    [temp,~] = min(tempf,[],'omitnan');
                    if any([isnan(temp),isempty(temp)])
                        continue;
                    end
                    numFreqPoints(k_pixels) = numel(tempf(~isnan(tempf)));
                    minFreqPoints(k_pixels) = min(tempf,[],'omitnan');
                    maxFreqPoints(k_pixels) = max(tempf,[],'omitnan');
                    if temp > minFreq
                        minFreq = temp;
                    end
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>lowfreqCutoff);
                    [temp,~] = max(tempf,[],'omitnan');
                    if isnan(temp)|| isempty(temp)
                        continue;
                    end
                    if temp < maxFreq
                        maxFreq = temp;
                    end
                end
                
                % Count orders of 10
                temp = minFreq;
                tempf = minFreq;
                magList = [];
                while temp < maxFreq
                    tempf = 10.^( floor( log10(temp) ) );
                    magList = horzcat(magList,tempf);
                    temp = temp + 200;
                end
                magList = unique(magList);
                
%                 % Cleverly make a list of numbers to sample at which
%                 % divide each order of magnitude into 10 samples.
%                 freqList = [];
%                 for k = 2:numel(magList)
%                     freqList = horzcat(freqList,linspace(magList(k-1),magList(k),10));
%                 end
%                 freqList = flip(unique(freqList)); % Process from high to low freq
                freqList = magList;
                
                fprintf('%s, # of Datapoints in Frequency: %d to %d; Median Number: %d, Frequency Range: %dHz to %dHz\n',...
                    fileLabels{j_dir},min(numFreqPoints,[],'omitnan'),...
                    max(numFreqPoints,[],'omitnan'),...
                    round(median(numFreqPoints,'all','omitnan')),...
                    min(minFreqPoints,[],'omitnan'),...
                    max(maxFreqPoints,[],'omitnan'));
                
                if plotKMeans
                    numObservations = NaN(size(resultsStruct.(varNames{j}).frequencyMap));
                    for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                        temp = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                        temp = temp(temp>=min(freqList) & temp<=max(freqList));
                        numObservations(k_pixels) = numel(temp);
                    end
                end
                
                if plotKMeans
                    
%                     % Cleverly make a list of numbers to sample at which
%                     % divide each order of magnitude into 10 samples.
%                     newFreqs = [];
%                     for k = 2:numel(freqList)
%                         newFreqs = horzcat(newFreqs,linspace(freqList(k-1),freqList(k),10));
%                     end
%                     newFreqs = unique(newFreqs);
                    newFreqs = magList;
                    
                    kmeansData = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),numel(newFreqs));
                    for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                        temp = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};
                        tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                        temp = temp(tempf>lowfreqCutoff);
                        tempf = tempf(tempf>lowfreqCutoff);
                        if isempty(temp)
                            continue;
                        end
                        
                        % Resample to known array of frequencies
%                         ids = ((tempf >= min(freqList)) & (tempf <= max(freqList)));
                        obsOut = interp1(tempf,temp,freqList,'makima',...
                            'extrap');
%                         try
%                             figure(hfig);
%                             clf;
%                         catch
%                             hfig = figure;
%                         end
%                         plot(newFreqs,obsOut,'r')
%                         hold on
%                         scatter(tempf,temp,'bo')
%                         hold off
                        
                        kmeansData(k_pixels,:) = obsOut;
                        
                    end                    
                    
                end
                                                        
                for k_freq = 1:numel(freqList)

                    % Make blank map data
                    mapDataStorage = NaN(mapSize);
                    mapDataLoss = NaN(mapSize);
                    mapDataAngle = NaN(mapSize);
                    mapDataRelaxance = NaN(mapSize);
                    mapDataError = NaN(mapSize);
                    mapDataTerms = NaN(mapSize);
                    mapDataHeight = NaN(mapSize);
                    mapDataKMeans = NaN(mapSize);
                    clf(mapPlotWindow)

                    % Position for the map
                    xc = 1;
                    yc = 0;
                    pixelLog = NaN(numel(mapDataStorage),2);

                    for k_pixels = 1:numel(mapDataStorage)
                        
                        % Get the current pixel position
                        if xc > mapSize(1)
                            xc = 1;
                            yc = yc + 1;
                        end
                        idx_pixel = sub2ind(mapSize,mapSize(2)-yc,xc);
                        pixelLog(k_pixels,:) = [mapSize(2)-yc,xc];
                        
%                         if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
%                             xc = xc + 1;
%                             continue;
%                         end

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
                            freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
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
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan'))
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataRelaxance(idx_pixel) = NaN;
                                mapDataError(idx_pixel) = NaN;
                                mapDataTerms(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                xc = xc + 1;
                            else
                                % Resample to known array of frequencies
%                                 ids = ((freq >= min(freqList)) & (freq <= max(freqList)));
                                mapDataStorage(idx_pixel) = interp1(freq,modelStorage,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataLoss(idx_pixel) = interp1(freq,modelLoss,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataAngle(idx_pixel) = interp1(freq,modelAngle,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataRelaxance(idx_pixel) = interp1(freq,modelRelaxance,freqList(k_freq),'makima',...
                                    'extrap');
                                                                
%                                 mapDataStorage(idx_pixel) = modelStorage(idx);
%                                 mapDataLoss(idx_pixel) = modelLoss(idx);
%                                 mapDataAngle(idx_pixel) = modelAngle(idx);
%                                 mapDataRelaxance(idx_pixel) = modelRelaxance(idx);
                                
                                mapDataError(idx_pixel) = modelErrorTime;
                                mapDataTerms(idx_pixel) = bestidx;
                                mapDataHeight(idx_pixel) = pixelHeight_cell{idx_pixel};
                                xc = xc + 1;
                            end
                            
                        else
                            
                            % We are just plotting the data captured
                            % DIRECTLY from the z-transform method.
                            freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                            modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                            modelAngle = atand(modelLoss./modelStorage);
                            modelRelaxance = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan'))
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataRelaxance(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                xc = xc + 1;
                            else
                                % Resample to known array of frequencies
%                                 ids = ((freq >= min(freqList)) & (freq <= max(freqList)));
                                mapDataStorage(idx_pixel) = interp1(freq,modelStorage,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataLoss(idx_pixel) = interp1(freq,modelLoss,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataAngle(idx_pixel) = interp1(freq,modelAngle,freqList(k_freq),'makima',...
                                    'extrap');
                                mapDataRelaxance(idx_pixel) = interp1(freq,modelRelaxance,freqList(k_freq),'makima',...
                                    'extrap');
                                
%                                 mapDataStorage(idx_pixel) = modelStorage(idx);
%                                 mapDataLoss(idx_pixel) = modelLoss(idx);
%                                 mapDataAngle(idx_pixel) = modelAngle(idx);

                                mapDataHeight(idx_pixel) = pixelHeight_cell{idx_pixel};
                                xc = xc + 1;
                            end
                        
                        end

                    end
                    
                    mapDataStorage(mapDataStorage == 0) = NaN;
                    mapDataLoss(mapDataLoss == 0) = NaN;
                    mapDataAngle(mapDataAngle == 0) = NaN;
                    
                    if plotKMeans
                        idxFreq = find(newFreqs == freqList(k_freq));
                        opts = statset('UseParallel',1);
                        idxK = kmeans(real(kmeansData(:,idxFreq)),kmeansBins,'Options',opts,'MaxIter',10000,...
                                'Display','off','Replicates',100);
                        for k_pixels = 1:numel(mapDataKMeans)
                            mapDataKMeans(pixelLog(k_pixels,1),pixelLog(k_pixels,2)) = idxK(k_pixels);
                        end
                        
                        figure(mapPlotKMeans)
                        colormap(mapColorName)
                        surf(X,Y,mapDataHeight,mapDataKMeans,'EdgeColor','interp')
                        hold on
                        title(sprintf('K-Means Clusters, %g Hz',freqList(k_freq)))
                        ylabel('Y Index')
                        xlabel('X Index')
                        xlim([1 mapSize(1)])
                        ylim([1 mapSize(2)])
                        cb = colorbar('Ticks',1:kmeansBins,...
                            'TickLabels',sprintfc('Bin %d',[1:kmeansBins]));
%                         caxis([1 kmeansBins]);
                        view(2)
                        hold off
                        
                        saveas(mapPlotKMeans,[path filesep fileLabels{j_dir} '-KMeansPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.fig'])
%                         saveas(mapPlotKMeans,[path filesep fileLabels{j_dir} '-KMeansPlot-' varNames{j} mapType '-' num2str(omegaList(k_omega)) 'Hz.jpg'])
                        print(mapPlotKMeans,[path filesep fileLabels{j_dir} '-KMeansPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.png'],'-dpng','-r300');

                    end
                    
                    if gradientOrder > 0
                        
                        for i_grad = 0:gradientOrder
                            
                            figure(mapPlotGradients)
                            colormap(mapColorName)
                            gradData = mapDataRelaxance;
                            GX = mapDataRelaxance;
                            GY = mapDataRelaxance;
                            ii = 0;
                            while ii < i_grad
                                [GX,GY] = gradient(abs(gradData));
                                gradData = sqrt(GX.^2+GY.^2);
    %                             gradData = GX+GY;
                                ii = ii + 1;
                            end
                            if i_grad == 0
                                surf(X,Y,mapDataHeight,abs(mapDataRelaxance),'EdgeColor','interp')
                            else
                                contourf(X,Y,abs(gradData),20);
                            end
                            hold on
                            if i_grad > 0
                                quiver(X,Y,GX,GY,4,'w-','linewidth',0.8)
                            end
                            title(sprintf('Relaxance Gradient (x%d), %g Hz',i_grad,freqList(k_freq)))
                            ylabel('Y Index')
                            xlabel('X Index')
                            xlim([1 mapSize(1)])
                            ylim([1 mapSize(2)])
                            cb = colorbar;
                            if i_grad > 1
                                cb.Ruler.TickLabelFormat=['%g' sprintf(' Pa/s^{%d}',i_grad)];
                            elseif i_grad == 1
                                cb.Ruler.TickLabelFormat=['%g Pa/s'];
                            else
                                cb.Ruler.TickLabelFormat=['%g Pa'];
                            end
                            view(2)
                            hold off

                            saveas(mapPlotGradients,[path filesep fileLabels{j_dir} '-MapGradient-Order' sprintf('%d',i_grad) '-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.fig'])
    %                         saveas(mapPlotGradients,[path filesep fileLabels{j_dir} '-MapGradient-Order' sprintf('%d',i_grad) '-' varNames{j} mapType '-' num2str(omegaList(k_omega)) 'Hz.jpg'])
                            print(mapPlotGradients,[path filesep fileLabels{j_dir} '-MapGradient-Order' sprintf('%d',i_grad) '-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.png'],'-dpng','-r300');
                        
                        end
                        
                    end
                    
                    figure(mapPlotWindow)
                    tiledlayout(1,3+plotModel, 'padding', 'none', 'TileSpacing', 'compact')
                    colormap(mapColorName)
                    nexttile
                    surf(X,Y,mapDataHeight,mapDataStorage,'EdgeColor','interp')
                    hold on
                    title(sprintf('Storage Modulus, %g Hz',freqList(k_freq)))
                    ylabel('Y Index')
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataStorage,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    cb.Ruler.TickLabelFormat='%g Pa';
                    caxis([0 climMax]);
                    view(2)
                    hold off

                    nexttile
                    surf(X,Y,mapDataHeight,mapDataLoss,'EdgeColor','interp')
                    hold on
                    title(sprintf('Loss Modulus, %g Hz',freqList(k_freq)))
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataLoss,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    cb.Ruler.TickLabelFormat='%g Pa';
                    caxis([0 climMax]);
                    view(2)
                    hold off
                    
                    nexttile
                    surf(X,Y,mapDataHeight,mapDataAngle,'EdgeColor','interp')
                    hold on
                    title(sprintf('Loss Angle, %g Hz',freqList(k_freq)))
                    xlabel('X Index')
                    xlim([1 mapSize(1)])
                    ylim([1 mapSize(2)])
%                     plotLims = prctile(mapDataAngle,plotRange,'all');
%                     zlim([plotLims(1)*0.8 plotLims(2)*1.2])
                    cb = colorbar;
                    cb.Ruler.TickLabelFormat='%g Deg';
                    caxis([0 90]);
                    view(2)
                    hold off

                    if plotModel
                        nexttile
                        surf(X,Y,mapDataHeight,mapDataTerms)
                        colormap(gca,'parula')
                        hold on
                        xlim([1 mapSize(1)])
                        ylim([1 mapSize(2)])
                        title(sprintf('Number of Terms'))
                        xlabel('X Index')
                        cb = colorbar;
                        set(cb,'YTick',1:numel(resultsStruct.(varNames{j}).bestParams))
                        view(2)
                        hold off
                    end

                    saveas(mapPlotWindow,[path filesep fileLabels{j_dir} '-MapPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.fig'])
                %     saveas(mapPlotWindow,[path filesep fileLabels{j_dir} '-MapPlot-' varNames{j} mapType '-' num2str(omegaList(k_omega)) 'Hz.jpg'])
                    print(mapPlotWindow,[path filesep fileLabels{j_dir} '-MapPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.png'],'-dpng','-r300');
                    
                end
                                    
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
            close all
            stillRunning = false;
        case 'Yes'
            close all
            clearvars -except stillRunning originalPath
    end
    
end
    
% Open the originally requested directory
winopen(originalPath);