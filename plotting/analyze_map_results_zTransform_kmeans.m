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
    dataType = 'storage'; % Type of kmeans clustering to perform
    figX = 10;
    figY = 10;
    nTicks = 5;
    plotModel = false;
    maxCol = 2;
    n_plots = 4+plotModel+plotIndentation;
    n_rows = ceil(n_plots/maxCol);
    n_cols = min([n_plots maxCol]);
    figWid = 375*n_cols;
    figHeight = 375*n_rows;
    mapColorName = 'turbo';
    climMax = 2e5; % Pa
    climHeight = 15e-6; % meters
    climInd = 2e-6; % meters
    stiffMax = 10*climMax; % Pa
    trimHeight = 50e-9;
    dFreq = 200; % Hz, step size between frames
    n_datapoints = 10;
    n_reps = 20; % number of kmeans replicates
    maxK = 10; % Max number of kmeans bins
    
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

    % Define error functions
    sse_global = @(data,model) sum((data-model).^2,'all');
    mse_global = @(data,model,n) sum((data-model).^2,'all')./(length(data)-n);
    
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
                
                minFreq = 0;
                maxFreq = Inf;
                for k_pixels = 1:numel(resultsStruct.(varNames{j}).frequencyMap)
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = min(tempf,[],'omitnan');
                    if any([isnan(temp),isempty(temp)]) || (numel(resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels})<n_datapoints)
                        continue;
                    end
                    if temp > minFreq
                        minFreq = temp;
                    end
                    tempf = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
                    tempf = tempf(tempf>0);
                    [temp,~] = max(tempf,[],'omitnan');
                    if isnan(temp)|| isempty(temp)
                        continue;
                    end
                    if temp < maxFreq
                        maxFreq = temp;
                    end
                end
                
                temp1 = min([minFreq,maxFreq]);
                temp2 = max([minFreq,maxFreq]);
                minFreq = 10;
                maxFreq = temp2;
                                
                % Count orders of 10
                temp = minFreq;
                tempf = minFreq;
                magList = [];
                while temp < maxFreq
                    tempf = 10.^( ( log10(temp) ) );
                    magList = horzcat(magList,tempf);
                    temp = temp + dFreq;
                end
                magList = unique(magList);
                
                pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
                pixelHeightArray = cell2mat(fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate));

                [minHeight,~] = min(pixelHeightArray);
                substrateCutoff = minHeight + trimHeight;
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel(pixelHeightArray);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list
                
                heightImg = zeros(mapSize);
                
                % Make blank map data
                mapDataHeight = NaN(mapSize);

                % Position for the map
                xc = 1;
                yc = 0;
                pixelLog = NaN(numel(mapDataHeight),2);
                
                clusteringData = NaN(numel(mapDataHeight),numel(magList));

                for k_pixels = 1:numel(mapDataHeight)

                    % Get the current pixel position
                    if xc > mapSize(1)
                        xc = 1;
                        yc = yc + 1;
                    end
                    idx_pixel = sub2ind(mapSize,mapSize(2)-yc,xc);
                    pixelLog(k_pixels,:) = [mapSize(2)-yc,xc];

                    heightImg(idx_pixel) = pixelHeightArray(idx_pixel);
                    
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

                    % We are just plotting the data captured
                    % DIRECTLY from the z-transform method.
                    modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                    modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
                    modelAngle = atand(modelLoss./modelStorage);
                        
                    % Resample to known array of frequencies
                    ids = ((freq >= min(magList)) & (freq <= max(magList)));
                    
                    if sum(ids,'all') < 2
                        xc = xc + 1;
                        continue;
                    end

%                     if hideSubstrate && (max([interp1(freq(ids),modelStorage(ids),magList,'makima',...
%                         NaN) interp1(freq(ids),modelLoss(ids),magList,'makima',...
%                         NaN) 0],[],'omitnan') > stiffMax)
% 
%                         xc = xc + 1;
%                         continue;
%                     end

                    kmeansInterp = [];

                    switch dataType
                        case 'force'
                            kmeansInterp = F_hz;
                        case 'indentation'
                            kmeansInterp = h_hz;
                        case 'storage'
                            kmeansInterp = modelStorage;
                        case 'loss'
                            kmeansInterp = modelLoss;
                        case 'angle'
                            kmeansInterp = modelAngle;
                    end

                    try
                        % Resample to known array of frequencies
                        obsOut = interp1(freq,kmeansInterp,magList,'makima',...
                            NaN);
                        clusteringData(k_pixels,:) = obsOut;
                    catch
                        % Do nothing
                    end
                    
                    mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                    xc = xc + 1;

                end

%                 opts = statset('UseParallel',1);
%                 tempfunc = @(x,k) kmeans(x,k,'Options',opts,'MaxIter',10000,...
%                         'Display','final','Replicates',n_reps);
                    
                opts = statset('UseParallel',1,...
                    'MaxIter',10000,...
                    'Display','final');
                tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
                    'Distance',@dtwf,...
                    'Replicates',n_reps);
                
                % Try all of the kmeans configurations
                eva = evalclusters(clusteringData,tempfunc,'CalinskiHarabasz',...
                    'klist',(1:maxK));
                idxK = tempfunc(clusteringData,eva.OptimalK);
                
                for k_pixels = 1:numel(mapDataHeight)
                    mapDataKMeans(pixelLog(k_pixels,1),pixelLog(k_pixels,2)) = idxK(k_pixels);
                end

                figure(mapPlotWindow)
                colormap(mapColorName)
                surf(X,Y,mapDataHeight,mapDataKMeans,'EdgeColor','interp')
                hold on
                title(sprintf('K-Means Clustering'))
                ylabel('Y Index')
                xlabel('X Index')
                xlim([1 mapSize(1)])
                ylim([1 mapSize(2)])
                cb = colorbar('Ticks',1:eva.OptimalK,...
                    'TickLabels',sprintfc('Bin %d',[1:eva.OptimalK]));
%                 caxis([1 eva.OptimalK]);
                view(2)
                hold off

                saveas(mapPlotWindow,[path filesep fileLabels{j_dir} '-KMeansPlot-' dataType '-' varNames{j} mapType '.fig'])
%                 saveas(mapPlotKMeans,[path filesep fileLabels{j_dir} '-KMeansPlot-' dataType '-' varNames{j} mapType '.jpg'])
                print(mapPlotWindow,[path filesep fileLabels{j_dir} '-KMeansPlot-' dataType '-' varNames{j} mapType '-' '.png'],'-dpng','-r300');

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