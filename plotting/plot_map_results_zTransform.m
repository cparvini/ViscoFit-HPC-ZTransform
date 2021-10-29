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
            clearvars timeErrorPlotNelder
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
        Files = dir([path '\*Results*zTransform*.mat']);
        
        if isempty(Files)
            error('The directory you selected does not contain a Z-Transform QI map. Please verify your FitResults file is in that directory and the filename contains "zTransform".');
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
                
                freqList = [];
                kmeansData = NaN(numel(resultsStruct.(varNames{j}).frequencyMap),numel(resultsStruct.(varNames{j}).frequencyMap{1}));
                for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                    temp = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    temp = 10.^(unique(floor(log10(temp))));
                    if any(isnan(temp)) || isempty(temp)
                        continue;
                    end
                    
                    % Remove frequencies to query where we have a curve
                    % that doesn't contain enough info (so all pixels can
                    % have data for our plots).
                    if k_pixels == 1
                        freqList = temp;
                    elseif numel(temp) < numel(freqList)
                        freqList = temp;
                    end
                end
                
                if plotKMeans
                    for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                        temp = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};
                        if isempty(temp)
                            continue;
                        end
                        kmeansData(k_pixels,:) = temp((temp >= min(freqList)) & (temp <= max(freqList)));
%                         kmeansData(:,k_pixels) = temp((temp >= min(freqList)) & (temp <= max(freqList)));
                    end
                    
                    opts = statset('UseParallel',1);
                    idx = kmeans(kmeansData,'Options',opts,'MaxIter',10000,...
                            'Display','final','Replicates',10);
                        
                    disp('pause')
                    
                end
                
                freqList = freqList(freqList>0);
                                                        
                for k_freq = 1:numel(freqList)

                    % Make blank map data
                    mapDataStorage = NaN(mapSize);
                    mapDataLoss = NaN(mapSize);
                    mapDataAngle = NaN(mapSize);
                    mapDataError = NaN(mapSize);
                    mapDataTerms = NaN(mapSize);
                    mapDataHeight = NaN(mapSize);
                    clf(mapPlotWindow)

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
                            
                            [~,idx] = min(abs(freq-freqList(k_freq)));
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan'))
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataError(idx_pixel) = NaN;
                                mapDataTerms(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                xc = xc + 1;
                                continue;
                            else
                                [~,idx] = min(abs(freq-freqList(k_freq)));
                                mapDataStorage(idx_pixel) = modelStorage(idx);
                                mapDataLoss(idx_pixel) = modelLoss(idx);
                                mapDataAngle(idx_pixel) = modelAngle(idx);
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
                            
                            
                            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan'))
                                mapDataStorage(idx_pixel) = NaN;
                                mapDataLoss(idx_pixel) = NaN;
                                mapDataAngle(idx_pixel) = NaN;
                                mapDataHeight(idx_pixel) = NaN;
                                xc = xc + 1;
                                continue;
                            else
                                [~,idx] = min(abs(freq-freqList(k_freq)));
                                mapDataStorage(idx_pixel) = modelStorage(idx);
                                mapDataLoss(idx_pixel) = modelLoss(idx);
                                mapDataAngle(idx_pixel) = modelAngle(idx);
                                mapDataHeight(idx_pixel) = pixelHeight_cell{idx_pixel};
                                xc = xc + 1;
                            end
                        
                        end

                    end
                    
%                     % Remove Outliers
%                     c = -1/(sqrt(2)*erfcinv(3/2));
%                     MAD_temp = 3.*c.*median(abs(mapDataStorage-median(mapDataStorage,'all')),'all');
%                     idx = abs(mapDataStorage-median(mapDataStorage,'all')) > MAD_temp;
% %                     MAD_temp = c.*median(abs(mapDataLoss-median(mapDataLoss,'all')),'all');
% %                     idx = idx + (abs(mapDataLoss-median(mapDataLoss,'all')) > MAD_temp);
%                     idx = logical(idx);
%                     mapDataStorage(idx) = NaN;
%                     mapDataLoss(idx) = NaN;
%                     mapDataAngle(idx) = NaN;
%                     mapDataError(idx) = NaN;
%                     mapDataTerms(idx) = NaN;

%                     plotRange = [10 90];
                    
                    mapDataStorage(mapDataStorage == 0) = NaN;
                    mapDataLoss(mapDataLoss == 0) = NaN;
                    mapDataAngle(mapDataAngle == 0) = NaN;

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

                    saveas(mapPlotWindow,[originalPath '\MapPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.fig'])
                %     saveas(mapPlotWindow,[originalPath '\mapPlot-' varNames{j} mapType '-' num2str(omegaList(k_omega)) 'Hz.jpg'])
                    print(mapPlotWindow,[originalPath '\MapPlot-' varNames{j} mapType '-' num2str(freqList(k_freq)) 'Hz.png'],'-dpng','-r300');

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
            stillRunning = false;
        case 'Yes'
            close all
            clearvars -except stillRunning originalPath
    end
    
end
    
% Open the originally requested directory
winopen(originalPath);