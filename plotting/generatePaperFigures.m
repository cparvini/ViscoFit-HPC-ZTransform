function [] = generatePaperFigures(originalPath,savePath,varargin)
%GENERATEPAPERFIGURES Create the Figures from the Main Manuscript
%   This function takes in a path argument and will subsequently create the
%   figures shown in the main manuscript. The datasets are exceedingly
%   large, and as such should be requested by email from the corresponding
%   author:
%
%   Alexander X. Cartagena-Rivera; alexander.cartagena-rivera@nih.gov
%
%   This script, and those contained within the repository, were written by
%   the first author:
%
%   Cameron H. Parvini; cameron.parvini@gmail.com
%   
%   While some of the details will have to be modified for users to
%   effectively run this code, they are relatively few. The specific cell
%   line names, dish numbers, and cell numbers would have to be changed in
%   particular, although the user could also provide "new" labels used in
%   their case to change what files are found using the search. The files
%   will, however, need to have already been processed using the other
%   functions in this directory so that the script can find the necessary
%   results datasets.

% User-Defined Settings
correctTilt = true;
hideSubstrate = true;
zeroSubstrate = true;
optimizeFlattening = true;
fillPixels = true;
showLabels = false;
showScaleBar = true;
clusterTarget = 'storage';
climMax = 2e5;
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
                        showLabels = varargin{i};                        
                    end
                case 7
                    if ~isempty(varargin{i})
                        showScaleBar = varargin{i};                        
                    end
                case 8
                    if ~isempty(varargin{i})
                        clusterTarget = varargin{i};                        
                    end
                case 9
                    if ~isempty(varargin{i})
                        if isa(varargin{i},'numeric')
                            climMax = varargin{i};
                        else
                            climMax = 2e5; % Pa
                            warning('You provided a non-numeric input for the stiffness plot limit to makeZAnimation(). Using the default value instead (2e5 Pa)');
                        end
                    end
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

%% Permanent Settings
errortype = 'sse';    
figX = 0;
figY = 0;
maxwid = get(0,'screensize');
maxheight = maxwid(4);
maxwid = maxwid(3);
figWid = min([maxwid maxheight]);
figHeight = figWid;
mapColorName = 'turbo';
boxSetting = 'on';
zMax = 15e-6; % meters, the JPK Nanowizard has a 15um piezo 
climInd = 1000e-9; % meters
stiffMax = 10*climMax; % Pa
trimHeight = 100e-9;
dFreq = 200; % Hz, step size between frames
n_steps = 100; % frames, number of frames per order of magnitude
n_datapoints = 10;
mediumFontSize = 14;
largeFontSize = 18;

%% Start Timer
tic

%% Figure 1 Data Selections
exampleCellType = {'HFF'};      % Cell Type to use
exampleCellDish = {1};          % Dish Number of Type
exampleCellNumber = {4};        % Cell Number in Dish
examplePixelNumber = {1000};    % Pixel to grab observables from
evalPt = 500;                   % Frequency for evaluating maps

%% Figure 1
%

fprintf('\nGenerating Figure 1...');
% Begin by finding the file we need for this figure and loading the data
% Search recursively (using "**") for our zTransform results files
dataFile = dir(fullfile(originalPath, '**',...
    ['*' exampleCellType{1} 'Dish' exampleCellDish{1} '*Cell' exampleCellNumber{1} '*Results*zTransform*.mat']));

if isempty(dataFile)
    error('The directory you selected does not contain the Z-Transform QI map indicated for Figure 1. Please verify your MapResults file is in that directory and the filename contains "zTransform".');
end

resultsStruct = load([dataFile.folder filesep dataFile.name],'-mat');
varNames = fields(resultsStruct);

% If there are any other analysis results...skip them! This is
% primarily for compatibility later, in case a new analysis type is
% proposed in the future.    
j = find(strcmpi(varNames,'zTransform'));
    
% Extract the size of the map from our results file
if isfield(resultsStruct.(varNames{j}),'mapSize')
    mapSize = resultsStruct.(varNames{j}).mapSize;
else
    mapSize = [128 128];
end

% In the case that we were able to extract the AFM-programmed map size,
% we will extract that here and use it for plotting our datasets on
% length-scaled axes. If this does not work, we could not determine
% how large the map was and we will use "X Index" and "Y Index"
% instead.
if isfield(resultsStruct.(varNames{j}),'scanSize')
    % We have the absolute map size!
    scanSize = resultsStruct.(varNames{j}).scanSize;
    xdataAbs = 0:(scanSize(1)/(mapSize(1)-1)):scanSize(1);
    ydataAbs = flip(0:(scanSize(2)/(mapSize(2)-1)):scanSize(2));
    [XA, YA] = meshgrid(xdataAbs,ydataAbs);
end

pixelHeight_cell = resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell;

% Axes meshgrid for scattering data. Here, we use the index instead of
% the absolute length/position.
xdata = 1:mapSize(1);
ydata = flip(1:mapSize(2));
[X, Y] = meshgrid(xdata,ydata);

% We'll need to calculate a consistent frequency array to use for
% comparing all of the pixels, but of course that requires taking a
% look at each pixel's frequency range FIRST so we can sample them all
% at the right points later.
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
            dFreq = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_steps);
        end

        if temp >= tempmax
            tempmax = temp*10;
            dFreq = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_steps);
        end

        tempf = 10.^( ( log10(temp) ) );
        magList = horzcat(magList,tempf);
        temp = temp + dFreq;
    end
end
magList = unique(magList);
freqList = flip(magList);

% Create an array representing each pixel's measured height, and
% correct that height according to the settings used during our
% analysis.
pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
if correctTilt
    temp = fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate,[],optimizeFlattening);
    pixelHeightArray = cell2mat(temp);
else
    pixelHeightArray = cell2mat(pixelHeight_cell);
end

% Remove pixels that are below our trim height to hopefully exclude the
% substrate. If not removed, the substrate will appear significantly
% stiffer than the cell and become saturated in the output plots.
[minHeight,~] = min(pixelHeightArray);
substrateCutoff = minHeight + trimHeight;
pixelsToRemove = false(size(pixelHeightArray));
pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

% Remove the pixels we want to keep from the list
pixelSkip = 1:numel(pixelHeightArray);
pixelSkip(~pixelsToRemove) = [];

% Starting position for the map
xc = 1;
yc = 0;

% First, if we don't have our indentation map stored (legacy files), we
% have to quickly calculate our observables to use during the analysis.
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

% Make blank map data
mapDataStorage = NaN(flip(mapSize));
mapDataLoss = NaN(flip(mapSize));
mapDataAngle = NaN(flip(mapSize));
mapDataRelaxance = NaN(flip(mapSize));
mapDataHeight = NaN(flip(mapSize));
mapDataInd = NaN(flip(mapSize));

pixelLog = NaN(numel(mapDataHeight),2);

for k_pixels = 1:numel(mapDataStorage)

    % Get the current pixel position
    if xc > mapSize(1)
        xc = 1;
        yc = yc + 1;
    end

    % Find our linear index for saving things later
    idx_pixel = sub2ind(flip(mapSize),mapSize(2)-yc,xc);
    pixelLog(k_pixels,:) = [mapSize(2)-yc,xc];
    
    % Skip the pixel if there is NaN data in the frequency vector
    if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
        xc = xc + 1;
        continue;
    end

    % Skip the pixel if we have decided to hide the substrate AND
    % this pixel is in our ignore list based on it's height
    if hideSubstrate && any(ismember(k_pixels,pixelSkip))
        xc = xc + 1;
        continue;
    end

    % If we had to calculate our observables, grab those values. If
    % not, great! We have the observables in our results file
    % directly.
    if ~isfield(resultsStruct.(varNames{j}),'indMap')
        F_hz = abs(F_hz_all{k_pixels});
        h_hz = abs(h_hz_all{k_pixels});
    else
        F_hz = abs(resultsStruct.(varNames{j}).forceMap{k_pixels});
        h_hz = abs(resultsStruct.(varNames{j}).indMap{k_pixels});
    end

    % Load and perform peak correction for the frequency vector.
    % This happens because of sampling not allowing us to truly
    % collect a "zero" frequency.
    freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
    [~,maxid] = max(F_hz);
    freqAdj = freq(maxid);
    freq = freq - freqAdj;

    % We will want to save our two-sided frequency vector for
    % later.
    if k_pixels == examplePixelNumber{1}
        freq_example = freq;
        F_hz_example = F_hz;
        h_hz_example = h_hz;
    end

    % We are just plotting the data captured
    % DIRECTLY from the z-transform method.
    modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
    modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
    modelAngle = atand(modelLoss./modelStorage);
    modelRelaxance = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};

    if (evalPt > max(freq,[],'omitnan')) || (evalPt < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq((freq >= min(freqList)) & (freq <= max(freqList)))) < 2

        mapDataStorage(idx_pixel) = NaN;
        mapDataLoss(idx_pixel) = NaN;
        mapDataAngle(idx_pixel) = NaN;
        mapDataRelaxance(idx_pixel) = NaN;
        mapDataHeight(idx_pixel) = NaN;
        mapDataInd(idx_pixel) = NaN;
        xc = xc + 1;

        % This is a problematic example, so we'll just give empty
        % data. If we see this in the output, we should re-run with
        % a different pixel example selected!
        if k_pixels == examplePixelNumber{1}
            h_t_example = NaN;
            F_t_example = NaN;
        end

    else

        ids = ((freq >= min(freqList)) & (freq <= max(freqList)));

        if hideSubstrate && (max([interp1(freq(ids),modelStorage(ids),evalPt,'makima',...
            NaN) interp1(freq(ids),modelLoss(ids),evalPt,'makima',...
            NaN) 0],[],'omitnan') > stiffMax)

            xc = xc + 1;
            continue;
        end

        mapDataStorage(idx_pixel) = interp1(freq(ids),modelStorage(ids),evalPt,'makima',...
            NaN);
        mapDataLoss(idx_pixel) = interp1(freq(ids),modelLoss(ids),evalPt,'makima',...
            NaN);
        mapDataAngle(idx_pixel) = interp1(freq(ids),modelAngle(ids),evalPt,'makima',...
            NaN);
        mapDataRelaxance(idx_pixel) = interp1(freq(ids),modelRelaxance(ids),evalPt,'makima',...
            NaN);

        % Has to be from the time domain
        % (normalization issues)
        t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
        h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
        mapDataInd(idx_pixel) = interp1(t_t,h_t,evalPt,'makima',...
            1e-12);
        mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);

        % We will want to save our two-sided frequency vector for
        % later.
        if k_pixels == examplePixelNumber{1}
            temph = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
            h_t_example = interp1(t_t,temph,evalPt,'makima',...
                1e-12);
            tempF = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
            F_t_example = interp1(t_t,tempF,evalPt,'makima',...
                1e-12);
        end

        xc = xc + 1;

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

end

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
else
    XPlot = X;
    xlab = 'X Index';
    xlims = [1 max(mapSize)];
    YPlot = Y;
    ylab = 'Y Index';
    ylims = [1 max(mapSize)];
end

% Quickly save these settings for Fig. 3
XPlotExample = XPlot;
xlabExample = xlab;
xlimsExample = xlims;
YPlotExample = YPlot;
ylabExample = ylab;
ylimsExample = ylims;

% Subfigure (a): Raw F-h Data and Z-Transform Example
% Clear old figures if they exist
if ~exist('mapPlotWindow','var')
    mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
else
    try
        figure(mapPlotWindow)
        clf
        mapPlotWindow.Position = [figX figY figWid figHeight];
    catch
        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
    end
end

% Plot raw example datasets
plot(h_t_example./1e-6,F_t_example./1e-9,'r-','LineWidth',5)
hold on
box(ax,boxSetting)
xlabel('Indentation [$$\\mum$$]','Interpreter','latex')
xlim([0 max(h_t_example)])
ylabel('Force [$$nN$$]','Interpreter','latex')
ylim([0 1.1*max(F_t_example)])
view(2)
pbaspect([1 1 1])
hold off

% Using exportgraphics() for Higher Quality
plotFile = [savePath filesep 'Fig1a-FDcurve'];
saveas(mapPlotWindow,[plotFile '.fig'])
exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
exportgraphics(ax,[plotFile '.png'],'Resolution',300);

% Clear old figures if they exist
if ~exist('mapPlotWindow','var')
    mapPlotWindow = figure('Position',[figX figY figWid*1.5 figHeight]);
else
    try
        figure(mapPlotWindow)
        clf
        mapPlotWindow.Position = [figX figY figWid*1.5 figHeight];
    catch
        mapPlotWindow = figure('Position',[figX figY figWid*1.5 figHeight]);
    end
end

% Plot Z-Transform example datasets
yyaxis left
plot(freq_example,F_hz_example,'r-','LineWidth',5)
hold on
box(ax,boxSetting)
xlabel('Frequency [$$Hz$$]','Interpreter','latex')
xlim([0 max(freq_example)])
ylabel('Force [$$Arb.$$]','Interpreter','latex')
ylim([0 1.1*max(F_hz_example)])
view(2)
hold off

yyaxis right
plot(freq_example,h_hz_example,'r-','LineWidth',5)
hold on
box(ax,boxSetting)
xlim([0 max(freq_example)])
ylabel('Indentation [$$Arb.$$]','Interpreter','latex')
ylim([0 1.1*max(h_hz_example)])
view(2)
hold off

pbaspect([1 1 1])

% Using exportgraphics() for Higher Quality
plotFile = [savePath filesep 'Fig1a-ZDistributions'];
saveas(mapPlotWindow,[plotFile '.fig'])
exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
exportgraphics(ax,[plotFile '.png'],'Resolution',300);

% Subfigure (b): Example Observables from QI Map
% Clear old figures if they exist
if ~exist('mapPlotWindow','var')
    mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
else
    try
        figure(mapPlotWindow)
        clf
        mapPlotWindow.Position = [figX figY figWid figHeight];
    catch
        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
    end
end

% Find the point that we used as our example in the other subfigure
xpos = pixelLog(examplePixelNumber{1},2); ypos = pixelLog(examplePixelNumber{1},1);

for i_plot = 1:5

    figure(mapPlotWindow)
    clf

    switch i_plot
        case 1
            mapData = mapDataHeight;
            plotTitle = 'Topography';
            saveLabel = 'Topography';
        
        case 2
            mapData = mapDataInd;
            plotTitle = 'Indentation';
            saveLabel = 'Indentation';
            
        case 3
            mapData = mapDataStorage;
            plotTitle = 'Storage Modulus';
            saveLabel = 'Storage';

        case 4
            mapData = mapDataLoss;
            plotTitle = 'Loss Modulus';
            saveLabel = 'Loss';

        case 5
            mapData = mapDataAngle;
            plotTitle = 'Loss Angle';
            saveLabel = 'Angle';

    end

    surf(XPlot,YPlot,mapDataHeight,mapData,'EdgeColor','interp')
    colormap(ax,mapColorName)
    hold on
    box(ax,boxSetting)
    title(plotTitle,'FontSize',largeFontSize)
    if showLabels
        xlabel(xlab)
        ylabel(ylab)
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end
    xlim(xlims)
    ylim(ylims)
    cb = colorbar;
    switch i_plot
        case 1
            caxis([0 zMax]);
            temp = (cb.Ticks' ./ 1e-6);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g \\mum',temp(ii));
            end
            scaleBarMax = zMax;
            
            % Add our marker for where the example plots came from!
            scatter3(xpos,ypos,mapDataHeight(xpos,ypos),25,'rs')
            
        case 2
            caxis([0 climInd]);
            temp = (cb.Ticks' ./ 1e-9);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g nm',temp(ii));
            end
            scaleBarMax = climInd;
        
        case 3
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
            end
            scaleBarMax = climMax;

        case 4
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
            end
            scaleBarMax = climMax;

        case 5
            cb.Ruler.TickLabelFormat='%g Deg';
            caxis([0 90]);
            scaleBarMax = 90;

    end
    
    if showScaleBar && exist('XA','var') && exist('YA','var')
        barsize = 10^ceil(log10(max(xlims))-1);
        xbar = [(1/20)*xlims(2) (1/20)*xlims(2)+barsize];
        ybar = ones(size(xbar))*(1/12.5)*xlims(2);
        zbar = ones(size(ybar))*scaleBarMax;
        plot3(xbar,ybar,zbar,'-k','LineWidth',5)
        text(xbar(1),ybar(1),zbar(1),sprintf('%d \\mum',barsize/1e-6),...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top',...
            'FontSize',mediumFontSize);
        
        xlim(xlims)
        ylim(xlims)
        zlim([0 scaleBarMax])
        view(2)
    end
        
    hold off
    
    view(2)
    pbaspect([1 1 1])
    drawnow

    % Create our filename
    plotFile = [savePath filesep 'Fig1b-' saveLabel];
    
    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(ax,[plotFile '.png'],'Resolution',300);

end

% Before moving on, while we are here, we are going to save some results
% for plotting later (Fig. 3).
mapDataHeightExample = mapDataHeight;
mapDataStorageExample = mapDataStorage;
cid = find(strcmp({resultsStruct.(varNames{j}).clusterData.clusterVar}, clusterTarget));
mapDataCluster2DExample = resultsStruct.(varNames{j}).clusterData(cid).clusterMap2D;
mapDataGlobalCluster2DExample = resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap2D;

fprintf('complete!\n');

%% Figure 2 Data Selections
exampleCellType = {'HEMa','HFF','A375P','A375M1','A375M2'};     % Cell Type to use
exampleCellDish = {1, 1, 1 ,1 ,1};                              % Dish Number of Type
exampleCellNumber = {1, 1, 1, 1, 1};                            % Cell Number in Dish
evalPt = [5e2, 1e3, 5e3, 10e3, 50e3];                           % Frequency for evaluating maps
dataOrder = {'ind','storage','loss','angle'};                   % Order (L > R, Top > Bottom) to plot observables in quad

%% Figure 2
% For this figure, we create compact grids of viscoelastic observables for
% each cell type at a variety of frequencies. There is quite a lot of data,
% so runtimes for this section may be large.

fprintf('\nGenerating Figure 2...');
for i_file = 1:numel(exampleCellType)

    % Begin by finding the file we need for this figure and loading the data
    % Search recursively (using "**") for our zTransform results files
    dataFile = dir(fullfile(originalPath, '**',...
        ['*' exampleCellType{i_file} 'Dish' exampleCellDish{i_file} '*Cell' exampleCellNumber{i_file} '*Results*zTransform*.mat']));

    if isempty(dataFile)
        error('The directory you selected does not contain the Z-Transform QI map indicated for Figure 1. Please verify your MapResults file is in that directory and the filename contains "zTransform".');
    end

    resultsStruct = load([dataFile.folder filesep dataFile.name],'-mat');
    varNames = fields(resultsStruct);

    % If there are any other analysis results...skip them! This is
    % primarily for compatibility later, in case a new analysis type is
    % proposed in the future.    
    j = find(strcmpi(varNames,'zTransform'));

    % Extract the size of the map from our results file
    if isfield(resultsStruct.(varNames{j}),'mapSize')
        mapSize = resultsStruct.(varNames{j}).mapSize;
    else
        mapSize = [128 128];
    end

    % In the case that we were able to extract the AFM-programmed map size,
    % we will extract that here and use it for plotting our datasets on
    % length-scaled axes. If this does not work, we could not determine
    % how large the map was and we will use "X Index" and "Y Index"
    % instead.
    if isfield(resultsStruct.(varNames{j}),'scanSize')
        % We have the absolute map size!
        scanSize = resultsStruct.(varNames{j}).scanSize;
        xdataAbs = 0:(scanSize(1)/(mapSize(1)-1)):scanSize(1);
        ydataAbs = flip(0:(scanSize(2)/(mapSize(2)-1)):scanSize(2));
        [XA, YA] = meshgrid(xdataAbs,ydataAbs);
    end

    pixelHeight_cell = resultsStruct.(varNames{j}).ViscoClass.pixelHeight_cell;

    % Axes meshgrid for scattering data. Here, we use the index instead of
    % the absolute length/position.
    xdata = 1:mapSize(1);
    ydata = flip(1:mapSize(2));
    [X, Y] = meshgrid(xdata,ydata);

    % Use our query list to loop through
    freqList = sort(evalPt);

    % Create an array representing each pixel's measured height, and
    % correct that height according to the settings used during our
    % analysis.
    pixelHeightArray = NaN(size([pixelHeight_cell{:}]));
    if correctTilt
        temp = fixMapTilt({mapSize},pixelHeight_cell,zeroSubstrate,[],optimizeFlattening);
        pixelHeightArray = cell2mat(temp);
    else
        pixelHeightArray = cell2mat(pixelHeight_cell);
    end

    % Remove pixels that are below our trim height to hopefully exclude the
    % substrate. If not removed, the substrate will appear significantly
    % stiffer than the cell and become saturated in the output plots.
    [minHeight,~] = min(pixelHeightArray);
    substrateCutoff = minHeight + trimHeight;
    pixelsToRemove = false(size(pixelHeightArray));
    pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

    % Remove the pixels we want to keep from the list
    pixelSkip = 1:numel(pixelHeightArray);
    pixelSkip(~pixelsToRemove) = [];

    % Starting position for the map
    xc = 1;
    yc = 0;

    % First, if we don't have our indentation map stored (legacy files), we
    % have to quickly calculate our observables to use during the analysis.
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
    
    for k_freq = 1:numel(freqList)
        
        % Make blank map data
        mapDataStorage = NaN(flip(mapSize));
        mapDataLoss = NaN(flip(mapSize));
        mapDataAngle = NaN(flip(mapSize));
        mapDataRelaxance = NaN(flip(mapSize));
        mapDataHeight = NaN(flip(mapSize));
        mapDataInd = NaN(flip(mapSize));

        pixelLog = NaN(numel(mapDataHeight),2);

        for k_pixels = 1:numel(mapDataStorage)

            % Get the current pixel position
            if xc > mapSize(1)
                xc = 1;
                yc = yc + 1;
            end

            % Find our linear index for saving things later
            idx_pixel = sub2ind(flip(mapSize),mapSize(2)-yc,xc);
            pixelLog(k_pixels,:) = [mapSize(2)-yc,xc];

            % Skip the pixel if there is NaN data in the frequency vector
            if any(isnan(resultsStruct.(varNames{j}).frequencyMap{k_pixels}))
                xc = xc + 1;
                continue;
            end

            % Skip the pixel if we have decided to hide the substrate AND
            % this pixel is in our ignore list based on it's height
            if hideSubstrate && any(ismember(k_pixels,pixelSkip))
                xc = xc + 1;
                continue;
            end

            % If we had to calculate our observables, grab those values. If
            % not, great! We have the observables in our results file
            % directly.
            if ~isfield(resultsStruct.(varNames{j}),'indMap')
                F_hz = abs(F_hz_all{k_pixels});
                h_hz = abs(h_hz_all{k_pixels});
            else
                F_hz = abs(resultsStruct.(varNames{j}).forceMap{k_pixels});
                h_hz = abs(resultsStruct.(varNames{j}).indMap{k_pixels});
            end

            % Load and perform peak correction for the frequency vector.
            % This happens because of sampling not allowing us to truly
            % collect a "zero" frequency.
            freq = resultsStruct.(varNames{j}).frequencyMap{k_pixels};
            [~,maxid] = max(F_hz);
            freqAdj = freq(maxid);
            freq = freq - freqAdj;

            % We are just plotting the data captured
            % DIRECTLY from the z-transform method.
            modelStorage = abs(real(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
            modelLoss = abs(imag(resultsStruct.(varNames{j}).relaxanceMap{k_pixels}));
            modelAngle = atand(modelLoss./modelStorage);
            modelRelaxance = resultsStruct.(varNames{j}).relaxanceMap{k_pixels};

            if (freqList(k_freq) > max(freq,[],'omitnan')) || (freqList(k_freq) < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq((freq >= min(freqList)) & (freq <= max(freqList)))) < 2

                mapDataStorage(idx_pixel) = NaN;
                mapDataLoss(idx_pixel) = NaN;
                mapDataAngle(idx_pixel) = NaN;
                mapDataRelaxance(idx_pixel) = NaN;
                mapDataHeight(idx_pixel) = NaN;
                mapDataInd(idx_pixel) = NaN;
                xc = xc + 1;

            else

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
                mapDataRelaxance(idx_pixel) = interp1(freq(ids),modelRelaxance(ids),freqList(k_freq),'makima',...
                    NaN);

                % Has to be from the time domain
                % (normalization issues)
                t_t = resultsStruct.(varNames{j}).ViscoClass.times_cell{k_pixels};
                h_t = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
                mapDataInd(idx_pixel) = interp1(t_t,h_t,freqList(k_freq),'makima',...
                    1e-12);
                mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);

                xc = xc + 1;

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

        end

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
        else
            XPlot = X;
            xlab = 'X Index';
            xlims = [1 max(mapSize)];
            YPlot = Y;
            ylab = 'Y Index';
            ylims = [1 max(mapSize)];
        end

        % Create our frequency quad plots
        % Clear old figures if they exist
        if ~exist('mapPlotWindow','var')
            mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
        else
            try
                figure(mapPlotWindow)
                clf
                mapPlotWindow.Position = [figX figY figWid figHeight];
            catch
                mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
            end
        end
        
        tiledlayout(2,2,...
            'padding', 'none', ...
            'TileSpacing', 'none')

        for i_plot = 1:4

            ax = nexttile;
            plotTarget = dataOrder{i_plot};

            switch plotTarget
                case 'ind'
                    mapData = mapDataInd;
                    plotTitle = 'Indentation';
                    saveLabel = 'Indentation';

                case 'storage'
                    mapData = mapDataStorage;
                    plotTitle = 'Storage Modulus';
                    saveLabel = 'Storage';

                case 'loss'
                    mapData = mapDataLoss;
                    plotTitle = 'Loss Modulus';
                    saveLabel = 'Loss';

                case 'angle'
                    mapData = mapDataAngle;
                    plotTitle = 'Loss Angle';
                    saveLabel = 'Angle';
                    
                case 'relaxance'
                    mapData = abs(mapDataRelaxance);
                    plotTitle = 'Relaxance';
                    saveLabel = 'Relaxance';

            end

            % We are hiding the colobars here.
            surf(XPlot,YPlot,mapDataHeight,mapData,'EdgeColor','interp')
            colormap(ax,mapColorName)
            hold on
            box(ax,boxSetting)
%             title(plotTitle,'FontSize',largeFontSize)
%             xlabel(xlab)
%             ylabel(ylab)
            set(gca,'YTickLabel',[],'XTickLabel',[])
            xlim(xlims)
            ylim(ylims)
%             cb = colorbar;
%             switch plotTarget
%                 case 'ind'
%                     caxis([0 climInd]);
%                     temp = (cb.Ticks' ./ 1e-9);
%                     for ii = 1:numel(temp)
%                        cb.TickLabels{ii} = sprintf('%g nm',temp(ii));
%                     end
% 
%                 case 'storage'
%                     caxis([0 climMax]);
%                     temp = (cb.Ticks' .* 1e-3);
%                     for ii = 1:numel(temp)
%                        cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
%                     end
% 
%                 case 'loss'
%                     caxis([0 climMax]);
%                     temp = (cb.Ticks' .* 1e-3);
%                     for ii = 1:numel(temp)
%                        cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
%                     end
% 
%                 case 'angle'
%                     cb.Ruler.TickLabelFormat='%g Deg';
%                     caxis([0 90]);
%                     
%                 case 'relaxance'
%                     caxis([0 climMax]);
%                     temp = (cb.Ticks' .* 1e-3);
%                     for ii = 1:numel(temp)
%                        cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
%                     end
% 
%             end

            hold off
            view(2)
            pbaspect([1 1 1])
            drawnow
            
        end
        
        % Create our filename
        plotFile = [savePath filesep 'Fig2-QuadPlot-' exampleCellType{i_file} '-' num2str(freqList(k_freq)) 'Hz'];
        
        % Using exportgraphics() for Higher Quality
        saveas(mapPlotWindow,[plotFile '.fig'])
        exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
        exportgraphics(ax,[plotFile '.png'],'Resolution',300);
    
    end
        
    % Plot our topography and save it
    % Clear old figures if they exist
    if ~exist('mapPlotWindow','var')
        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
    else
        try
            figure(mapPlotWindow)
            clf
            mapPlotWindow.Position = [figX figY figWid figHeight];
        catch
            mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
        end
    end

    plotTitle = 'Topography';
    saveLabel = 'Topography';

    surf(XPlot,YPlot,mapDataHeight,mapDataHeight,'EdgeColor','interp')
    colormap(ax,mapColorName)
    hold on
    box(ax,boxSetting)
    title(plotTitle,'FontSize',largeFontSize)
    if showLabels
        xlabel(xlab)
        ylabel(ylab)
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end
    xlim(xlims)
    ylim(ylims)
    cb = colorbar;
    caxis([0 zMax]);
    temp = (cb.Ticks' ./ 1e-6);
    for ii = 1:numel(temp)
       cb.TickLabels{ii} = sprintf('%g \\mum',temp(ii));
    end
    hold off

    view(2)
    pbaspect([1 1 1])
    drawnow

    % Create our filename
    plotFile = [savePath filesep 'Fig2-' saveLabel '-' exampleCellType{i_file}];

    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(ax,[plotFile '.png'],'Resolution',300);
    
    clearvars resultsStruct

end

fprintf('complete!\n');

%% Figure 3
% In this figure, we show how applying an unsupervised machine learning
% algorithm (kmedoids clustering) allows segmentation according to
% meaningful nanoscale viscoelastic properties near the surface of cells.
% The example from Figure 1 is now reproduced to show how the clustering
% algorithm has segmented the results according to Storage Modulus when
% looking at the individual maps (isolated clustering) and the entire
% population of cell type (global clustering).

% This visualizes the results we will summarize in Figure 4 at a larger
% level. Thankfully, we have saved our global results right next to our
% isolated results so accessing them is simple.

% Here we will make four panels in individual plots:
% 1. Surface Topography of the example cell
% 2. Storage Modulus at a specific frequency of the example cell
% 3. The results (bins) from isolated clustering of this map
% 4. The results (bins) from global clustering of all example cell's type
% Subfigure (b): Example Observables from QI Map
% Clear old figures if they exist

fprintf('\nGenerating Figure 3...');

if ~exist('mapPlotWindow','var')
    mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
else
    try
        figure(mapPlotWindow)
        clf
        mapPlotWindow.Position = [figX figY figWid figHeight];
    catch
        mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
    end
end

% Loop through our sublots to create them
for i_plot = 1:5

    figure(mapPlotWindow)
    clf

    switch i_plot
        case 1
            mapData = mapDataHeightExample;
            plotTitle = 'Topography';
            saveLabel = 'Topography';
            
        case 2
            mapData = mapDataStorageExample;
            plotTitle = 'Storage Modulus';
            saveLabel = 'Storage';

        case 3
            mapData = mapDataCluster2DExample;
            plotTitle = 'Locally-Observed Clusters';
            saveLabel = 'SoloCluster';

        case 4
            mapData = mapDataGlobalCluster2DExample;
            plotTitle = 'Globally-Observed Clusters';
            saveLabel = 'GlobalCluster';

    end

    surf(XPlotExample,YPlotExample,mapDataHeightExample,...
        mapData,'EdgeColor','interp')
    colormap(ax,mapColorName)
    hold on
    box(ax,boxSetting)
    title(plotTitle,'FontSize',largeFontSize)
    if showLabels
        xlabel(xlabExample)
        ylabel(ylabExample)
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end
    xlim(xlimsExample)
    ylim(ylimsExample)
    cb = colorbar;
    switch i_plot
        case 1
            caxis([0 zMax]);
            temp = (cb.Ticks' ./ 1e-6);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g \\mum',temp(ii));
            end
            
            % Add our marker for where the example plots came from!
            scatter3(xpos,ypos,mapDataHeight(xpos,ypos),25,'rs')
            
        case 2
            caxis([0 climInd]);
            temp = (cb.Ticks' ./ 1e-9);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g nm',temp(ii));
            end
        
        case 3
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
            end

        case 4
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d kPa',temp(ii));
            end

        case 5
            cb.Ruler.TickLabelFormat='%g Deg';
            caxis([0 90]);

    end
        
    hold off
    
    view(2)
    pbaspect([1 1 1])
    drawnow

    % Create our filename
    plotFile = [savePath filesep 'Fig3-' saveLabel];
    
    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(ax,[plotFile '.png'],'Resolution',300);

end

fprintf('complete!\n');

%% Figure 4 Data Selections
clusterCellTypes = {'HEMa','HFF','A375P','A375M1','A375M2'};    % Cell Type to use
evalPt = 5e2;                                                   % Frequency for evaluating maps
colorList = hsv(10);                                            % Colors for Plotting
nBins = 100;                                                    % Number of bins for the histogram

%% Figure 4
% This figure shows a histogram comparison of the cell types by several
% details:
% 1. Cluster (according to color)
% 2. Cell Type (according to subplot)
% 3. Storage Modulus Value at a specific frequency (X-axis Value)
% 4. Number of pixels in that bin (Y-Axis Value)

% There are a lot of features in this plot, so choosing how to visualize
% them all is tricky. However, this allows readers to compare cell types
% (each plot) according to similar colors---e.g. they can look at the red
% cluster for two different cell types and see how they differ. This will
% also give a clear and obvious way to differentiate between types that
% have different numbers of bins. To do so, we find the medoid of each bin
% and then sort them according to their value. Then, we assign a color for
% each bin such that the lowest bin colors match, as do the second lowest
% and so-on for each cell type. This plot will be a stack of five cell
% types, each with the same X-axis limits (0 to 200 kPa, as with our other
% plots).

% Because we lack the dimensions to visualize different cells on top of one
% another, we will use the GLOBAL clustering results and count the number
% of pixels inside each bin at a cell-type level.

fprintf('\nGenerating Figure 4...');

for i_dir = 1:length(clusterCellTypes)
    
    % Search recursively (using "**") for our zTransform results files
    Files = dir(fullfile(originalPath, '**',...
        ['*' clusterCellTypes{i_dir} '*Results*zTransform*.mat']));

    if isempty(Files)
        error('The directory you selected and its subdirectories do not contain a Z-Transform QI map for cell type No. %d. Please verify your results file is in or under that directory and the filename contains "zTransform".',i_dir);
    end

    % Since we are analyzing many maps, we must create our frequency
    % array ahead of time.
    minFreq = Inf;
    maxFreq = 0;
    mapsizes = NaN(numel(Files),1);

    for j_dir = 1:length(Files)
        
        vars = whos('-file',[Files(j_dir).folder filesep Files(j_dir).name]);
        varNames = {vars.name};
        
        % If there are any other analysis results...skip them! This is
        % primarily for compatibility later, in case a new analysis type is
        % proposed in the future.    
        j = find(strcmpi(varNames,'zTransform'));
        
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
    timeList = 1./freqList;
                
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
%     clusteringDataGlobal = NaN(bi(end),numel(magList));
    XDataGlobal = NaN(bi(end),1);
    YDataGlobal = NaN(bi(end),1);
    pixelLogGlobal = NaN(bi(end),2);
    mapIDglobal = NaN(bi(end),1);
                
    % Now, loop through the files and create our dataset
    fprintf('\nPopulating our global clustering dataset No. %d of %d...',i_dir,length(clusterCellTypes));
    for j_dir = 1:length(Files)

        resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
        varNames = fields(resultsStruct);

        % If there are any other analysis results...skip them! This is
        % primarily for compatibility later, in case a new analysis type is
        % proposed in the future.    
        j = find(strcmpi(varNames,'zTransform'));
        
        if isfield(resultsStruct.(varNames{j}),'mapSize')
            mapSize = resultsStruct.(varNames{j}).mapSize;
        else
            mapSize = [128 128];
        end

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

%         clusteringData = NaN(numel(pixelHeightArray),numel(magList));
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
%                 obsOut = interp1(xinterp,clusterInterp,obsList,'makima',...
%                     NaN);
%                 clusteringData(k_pixels,:) = obsOut;
                obsOut = interp1(xinterp,clusterInterp,evalPt,'makima',...
                    NaN);
                histXData(k_pixels) = obsOut;
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

        cid = find(strcmp({resultsStruct.(varNames{j}).clusterData.clusterVar}, clusterTarget));
        histYData = cell2mat(resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap)';

        % We want to order the clusters for this map according to
        % increasing magnitude (softest to largest) so all of our
        % clusters are aligned between the maps.
        uniqueBins = unique(histYData);
        medVal = NaN(size(uniqueBins));
        for k_align = 1:length(uniqueBins)

            % Grab our bin values at the frequency we will be plotting!
            % Then, we can find the median of that bin and use it to
            % sort our bin labels in increasing order.
            tempData = histXData(histYData == uniqueBins(k_align));
            medVal(k_align) = median(tempData,'omitnan');

        end
        [~,idold] = sort(medVal);
        idnew = uniqueBins(idold);
        histYDataTemp = NaN(size(histYData));
        for k_align = 1:length(idold)

            % Reorganize our data according to the new order
            histYDataTemp(histYData == idold(k_align)) = idnew(k_align);

        end
        histYData = histYDataTemp;
            
        % Concatenate our Clustering Data from this map and save
        % the indices for this particular file
        pixelLogGlobal(ai(j_dir):bi(j_dir),:) = pixelLog;
        mapIDglobal(ai(j_dir):bi(j_dir)) = j_dir;
%         clusteringDataGlobal(ai(j_dir):bi(j_dir),:) = clusteringData;
        XDataGlobal(ai(j_dir):bi(j_dir),:) = histXData;
        YDataGlobal(ai(j_dir):bi(j_dir),:) = histYData;
            
        clearvars resultsStruct mapSize histXData histYData histYDataTemp ...
            pixelLog

    end
        
    fprintf('complete!\n');
        
    % Now, go through and save the results to each output file and make
    % our plots!
    fprintf('\nGenerating our histogram plot...');
    
    % Clear old figures if they exist
    if ~exist('mapPlotWindow','var')
        mapPlotWindow = figure('Position',[figX figY figWid*2 figHeight]);
    else
        try
            figure(mapPlotWindow)
            clf
            mapPlotWindow.Position = [figX figY figWid*2 figHeight];
        catch
            mapPlotWindow = figure('Position',[figX figY figWid*2 figHeight]);
        end
    end
    
    % Now, begin making our histogram
    uniqueBins = unique(YDataGlobal);
    hi = gobjects(numel(uniqueBins),1);
    binEdges = linspace(0,climMax,nBins+1);
    hold on
    for i_plot = 1:numel(uniqueBins)
        % Create our histogram
        hi(i_plot) = histogram(histXData(histYData == uniqueBins(i_plot)),...
            'BinEdges', binEdges,...
            'FaceColor', colorList(i_plot,:),...
            'FaceAlpha', 0.5,...
            'EdgeColor', [0 0 0],...
            'EdgeAlpha', 0.5);
        
        % Plot a normal distribution over the top
        [mu,sigma] = normfit(hi(i_plot).Values);
        xvals = binEdges(1:end-1) + diff(binEdges)/2;
        yvals = exp(-(hi(i_plot).Values-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
        plot(xvals,yvals,'Color',colorList(i_plot,:),'LineWidth',5)
    end
    grid on
    xlim([0 climMax])
    
    if i_dir == length(clusterCellTypes)
        switch clusterTarget
            case 'force'
                xlab = sprintf('Force [N]');
            case 'indentation'
                xlab = sprintf('Indentation [\\mum]');
            case 'storage'
                xlab = sprintf('Storage Modulus [Pa]');
            case 'loss'
                xlab = sprintf('Loss Modulus [Pa]');
            case 'angle'
                xlab = sprintf('Loss Angle [deg]');
            case 'relaxance'
                xlab = sprintf('Relaxance [Pa]');
        end
        xlabel(xlab, 'FontSize', mediumFontSize)
        ylabel('Number of Pixels', 'FontSize', mediumFontSize)
    else
        ylabel('Number of Pixels', 'FontSize', mediumFontSize)
    end
    
    hold off
    
    % Create our filename
    plotFile = [savePath filesep 'Fig4-' clusterCellTypes{i_dir}];
    
    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(ax,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(ax,[plotFile '.png'],'Resolution',300);
    
    fprintf('complete for %s!\n',clusterCellTypes{i_dir});
    
end

fprintf('Figure 4 complete!\n');

%% Stop Timer
totalTime = toc;
totalHr = floor(totalTime/3600);
totalMin = floor((totalTime - totalHr*3600)/60);
totalSec = floor(totalTime - totalHr*3600 - totalMin*60);

%% Open output directory
fprintf('\nManuscript Figure Generation Complete!\n');
fprintf('\nTotal Time (approx.): %d hours, %d minutes %d seconds\n',totalHr,totalMin,totalSec);
% fprintf('\nTotal Time (approx.) was %d:%02d:%02d\n',totalHr,totalMin,totalSec);
winopen(savePath)

end