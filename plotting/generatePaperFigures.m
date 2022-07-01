function [] = generatePaperFigures(clusterPath,highResPath,savePath,varargin)
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
hideSubstrate = true;
fillPixels = true;
logSteps = true;
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
                        showLabels = varargin{i};                        
                    end
                case 5
                    if ~isempty(varargin{i})
                        showScaleBar = varargin{i};                        
                    end
                case 6
                    if ~isempty(varargin{i})
                        clusterTarget = varargin{i};                        
                    end
                case 7
                    if ~isempty(varargin{i})
                        if isa(varargin{i},'numeric')
                            climMax = varargin{i};
                        else
                            climMax = 2e5; % Pa
                            warning('You provided a non-numeric input for the stiffness plot limit to makeZAnimation(). Using the default value instead (2e5 Pa)');
                        end
                    end
                otherwise
                    fprintf('Passed additional parameters to generatePaperFigures() which were not used.');
            end
        end
    end
end

%% Permanent Settings
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
climInd = 200e-9; % meters
climForce = 1e-9; % nN
stiffMax = 10*climMax; % Pa
trimHeight = 100e-9;
dFreq = 200; % Hz, step size between frames
n_steps = 100; % frames, number of frames per order of magnitude
n_datapoints = 10;
tinyFontSize = 16;
smallFontSize = 28;
mediumFontSize = 36;
largeFontSize = 48;

% Placeholders if our file doesn't have these settings
correctTilt = true;
zeroSubstrate = true;
optimizeFlattening = true;

%% Start Timer
tic

%% Figure 1 Data Selections
% Note that we will be using a high resolution map for visualizations in
% Figure 2. For figures 1 and 3, we need to use a low resolution map
% because we have excluded the high-res maps from our global analysis and
% we want to show the isolated vs. global in figure 3.

% Make sure that the cell type, dish number, and cell number 
% actually represent a valid selection for this subgroup.

exampleCellType = {'HFF'};          % Cell Type to use
exampleCellDish = {1};              % Dish Number of Type
exampleCellNumber = {4};            % Cell Number in Dish
examplePixelNumber = {8000};        % Pixel to grab observables from
evalPt = 1000;                      % Frequency for evaluating maps

%% Figure 1
fprintf('Generating Figure 1...');
% Begin by finding the file we need for this figure and loading the data
% Search recursively (using "**") for our zTransform results files
dataFile = dir(fullfile(clusterPath, '**',...
    ['*' num2str(exampleCellType{1}) 'Dish' num2str(exampleCellDish{1})...
    '*Cell' num2str(exampleCellNumber{1}) '_*Results*zTransform*.mat']));

if isempty(dataFile)
    error('The directory you selected does not contain the Z-Transform QI map indicated for Figure 1. Please verify your MapResults file is in that directory and the filename contains "zTransform".');
end

resultsStruct = load([dataFile.folder filesep dataFile.name],'-mat');
varNames = fields(resultsStruct);

% If there are any other analysis results...skip them! This is
% primarily for compatibility later, in case a new analysis type is
% proposed in the future.
j = find(contains(varNames,'zTransform'));
    
% Extract the size of the map from our results file
if isfield(resultsStruct.(varNames{j}),'mapSize')
    mapSize = resultsStruct.(varNames{j}).mapSize;
else
    mapSize = [128 128];
end

% Grab some relevant settings
try
    correctTilt = resultsStruct.(varNames{j}).correctTilt;
    zeroSubstrate = resultsStruct.(varNames{j}).zeroSubstrate;
    optimizeFlattening = resultsStruct.(varNames{j}).optimizeFlattening;
catch
    warning('The variables correctTilt, zeroSubstrate, and optimizeFlattening were not found in the results file. Using hard-coded settings in generatePaperFigures() (all set to true).');
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
[minHeight,~] = min(pixelHeightArray(pixelHeightArray>0));
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
mapDataForce = NaN(flip(mapSize));
heightImg = zeros(flip(mapSize));

pixelLog = NaN(numel(mapDataHeight),2);
heightImg = rot90(reshape(pixelHeightArray,mapSize),1);

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
    
    % We will want to save our two-sided frequency vector for
    % later.
    if k_pixels == examplePixelNumber{1}
        freq_example = freq;
        F_hz_example = F_hz;
        h_hz_example = h_hz;
        Estorage_hz_example = modelStorage;
        Eloss_hz_example = modelLoss;
    end
    
    if (evalPt > max(freq,[],'omitnan')) || (evalPt < min(freq,[],'omitnan')) || any(isnan(freq)) || numel(freq((freq >= min(freqList)) & (freq <= max(freqList)))) < 2

        mapDataStorage(idx_pixel) = NaN;
        mapDataLoss(idx_pixel) = NaN;
        mapDataAngle(idx_pixel) = NaN;
        mapDataRelaxance(idx_pixel) = NaN;
        mapDataHeight(idx_pixel) = NaN;
        mapDataInd(idx_pixel) = NaN;
        mapDataForce(idx_pixel) = NaN;
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
        F_t = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
        mapDataInd(idx_pixel) = interp1(t_t,h_t,1./evalPt,'makima',...
            1e-12);
        mapDataForce(idx_pixel) = interp1(t_t,F_t,1./evalPt,'makima',...
            1e-12);
%         mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
        mapDataHeight(idx_pixel) = heightImg(mapSize(2)-yc,xc);
        

        % We will want to save our two-sided frequency vector for
        % later.
        if k_pixels == examplePixelNumber{1}
            temph = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{k_pixels};
%             h_t_example = interp1(t_t,temph,1./freq(ids),'makima',...
%                 1e-12);
            h_t_example = temph;
            tempF = resultsStruct.(varNames{j}).ViscoClass.forces_cell{k_pixels};
%             F_t_example = interp1(t_t,tempF,1./freq(ids),'makima',...
%                 1e-12);
            F_t_example = tempF;
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
mapDataForce(mapDataForce==0) = NaN;

if exist('XA','var') && exist('YA','var')
    XPlot = XA./1e-6;
    xlab = sprintf('X Position [\\mum]');
    xlims = [0 max(scanSize)./1e-6];
    XAExample = XPlot;
    YPlot = YA./1e-6;
    ylab = sprintf('Y Position [\\mum]');
    ylims = [0 max(scanSize)./1e-6];
    YAExample = YPlot;
    scanSizeExample = scanSize;
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
ax = gca;
plot(h_t_example,F_t_example,'r-','LineWidth',5)
% plot(h_t_example,smooth(F_t_example,floor(numel(F_t_example)*0.1)),'r-','LineWidth',5)
hold on
grid on
box(ax,boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
set( findall(gcf, '-property', 'fontsize'), 'fontsize', smallFontSize)
xlabel('Indentation [$$nm$$]','Interpreter','latex','FontSize',mediumFontSize)
xlim([0 max(h_t_example)])
ylabel('Force [$$pN$$]','Interpreter','latex','FontSize',mediumFontSize)
ylim([0 1.1*max(F_t_example)])
view(2)

temp = (xticks' ./ 1e-9);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
xticklabels(TickLabels)

temp = (yticks' ./ 1e-12);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
yticklabels(TickLabels)

hold off
pbaspect([1 1 1])

% Using exportgraphics() for Higher Quality
plotFile = [savePath filesep 'Fig1a-FDcurve'];
saveas(mapPlotWindow,[plotFile '.fig'])
exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

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
ax = gca;
plot(freq_example,F_hz_example,'LineWidth',8)
% plot(freq_example,smooth(F_hz_example,floor(numel(F_hz_example)*0.1)),'LineWidth',8)
hold on
grid on
box(ax,boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
set( findall(gca, '-property', 'fontsize'), 'fontsize', smallFontSize)
xlabel('Frequency [$$kHz$$]','Interpreter','latex','FontSize',mediumFontSize)
xlim([0 max(freq_example)])
ylabel('Force [$$Arb.$$]','Interpreter','latex','FontSize',mediumFontSize)
% ylim([0 1.1*max(F_hz_example)])
ylim([0 1.1*max(smooth(F_hz_example,floor(numel(F_hz_example)*0.1)))])
view(2)
hold off

yyaxis right
ax = gca;
plot(freq_example,h_hz_example,'LineWidth',5)
% plot(freq_example,smooth(h_hz_example,floor(numel(h_hz_example)*0.1)),'LineWidth',5)
hold on
grid on
box(ax,boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
xlim([0 max(freq_example)])
ylabel('Indentation [$$Arb.$$]','Interpreter','latex','FontSize',mediumFontSize)
ylim([0 1.1*max(h_hz_example)])
ylim([0 1.1*max(smooth(h_hz_example,floor(numel(h_hz_example)*0.1)))])
view(2)
hold off

temp = (xticks' ./ 1e3);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
xticklabels(TickLabels)

pbaspect([1 1 1])

% Using exportgraphics() for Higher Quality
plotFile = [savePath filesep 'Fig1a-ZDistributions'];
saveas(mapPlotWindow,[plotFile '.fig'])
exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

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

% Plot Z-Transform viscoelastic observables example datasets
% First, find ylim max because both axes have the same units
tempMax = [1.1*max(Estorage_hz_example) 1.1*max(Eloss_hz_example)];
pbaspect([1 1 1])

% Now, make the plot
yyaxis left
ax = gca;
plot(freq_example,Estorage_hz_example,'LineWidth',8)
% plot(freq_example,smooth(Estorage_hz_example,floor(numel(Estorage_hz_example)*0.1)),'LineWidth',8)
hold on
grid on
box(ax,boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
set( findall(gca, '-property', 'fontsize'), 'fontsize', smallFontSize)
xlabel('Frequency [$$kHz$$]','Interpreter','latex','FontSize',mediumFontSize)
xlim([0 max(freq_example)])
ylabel('Storage Modulus [$$kPa$$]','Interpreter','latex','FontSize',mediumFontSize)
ylim([0 max(tempMax)])

temp = (yticks' ./ 1e3);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
yticklabels(TickLabels)

view(2)
hold off

yyaxis right
ax = gca;
plot(freq_example,Eloss_hz_example,'LineWidth',8)
% plot(freq_example,smooth(Eloss_hz_example,floor(numel(Eloss_hz_example)*0.1)),'LineWidth',5)
hold on
grid on
box(ax,boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',2)

xlim([0 max(freq_example)])
ylabel('Loss Modulus [$$kPa$$]','Interpreter','latex','FontSize',mediumFontSize)
ylim([0 max(tempMax)])
view(2)
hold off

temp = (yticks' ./ 1e3);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
yticklabels(TickLabels)

temp = (xticks' ./ 1e3);
TickLabels = cell(size(temp));
for ii = 1:numel(temp)
   TickLabels{ii} = sprintf('%d',round(temp(ii)));
end
xticklabels(TickLabels)

% Using exportgraphics() for Higher Quality
plotFile = [savePath filesep 'Fig1a-ZViscoModuli'];
saveas(mapPlotWindow,[plotFile '.fig'])
exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

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
    ax = gca;
    colormap(ax,mapColorName)
    hold on
    box(ax,boxSetting)
    set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
    set( findall(gca, '-property', 'fontsize'), 'fontsize', smallFontSize)
    if showLabels
        title(plotTitle,'FontSize',largeFontSize)
        xlabel(xlab,'FontSize',mediumFontSize)
        ylabel(ylab,'FontSize',mediumFontSize)
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end
    xlim(xlims)
    ylim(ylims)
%     axis('equal')
    cb = colorbar;
    switch i_plot
        case 1
            caxis([0 zMax]);
            temp = (cb.Ticks' ./ 1e-6);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g',temp(ii)); % \\mum
            end
            scaleBarMax = zMax;
            
            % Add our marker for where the example plots came from!
            scatter3(XPlot(ypos,xpos),YPlot(ypos,xpos),mapDataHeight(ypos,xpos)*1.1,50,'rs','filled')
            
        case 2
            caxis([0 climInd]);
            temp = (cb.Ticks' ./ 1e-9);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g',temp(ii)); % nm
            end
            scaleBarMax = climInd;
        
        case 3
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d',round(temp(ii))); % kPa
            end
            scaleBarMax = climMax;

        case 4
            caxis([0 climMax]);
            temp = (cb.Ticks' .* 1e-3);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d',round(temp(ii))); % kPa
            end
            scaleBarMax = climMax;

        case 5
            cb.Ruler.TickLabelFormat='%g'; % Deg
            caxis([0 90]);
            cb.Ticks = [0 30 60 90];
            scaleBarMax = 90;

    end
    cb.FontSize = smallFontSize;
    
    if showScaleBar && exist('XA','var') && exist('YA','var')
        barsize = 10^ceil(log10(max(xlims))-1);
        xbar = [(1/20)*xlims(2) (1/20)*xlims(2)+barsize];
        ybar = ones(size(xbar))*(1/12.5)*xlims(2);
        zbar = ones(size(ybar))*scaleBarMax;
        plot3(xbar,ybar,zbar,'-k','LineWidth',5)
        text(xbar(1),ybar(1),zbar(1),sprintf('%d \\mum',barsize),...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top',...
            'FontSize',smallFontSize);
        
        xlim(xlims)
        ylim(xlims)
%         zlim([0 scaleBarMax])
        view(2)
    end
       
    if any(diff(yticks) ~= barsize)
        temp = yticks;
        yticks(ax,temp(1):barsize:temp(end))
    end
    if any(diff(xticks) ~= barsize)
        temp = xticks;
        xticks(ax,temp(1):barsize:temp(end))
    end
    
    hold off
    
    view(2)
    pbaspect([1 1 1])
    drawnow

    % Create our filename
    plotFile = [savePath filesep 'Fig1b-' saveLabel];
    
    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

end

% Before moving on, while we are here, we are going to save some results
% for plotting later (Fig. 3).
mapDataHeightExample = mapDataHeight;
switch clusterTarget
    case 'force'
        mapDataTargetExample = mapDataForce;
    case 'indentation'
        mapDataTargetExample = mapDataInd;
    case 'storage'
        mapDataTargetExample = mapDataStorage;
    case 'loss'
        mapDataTargetExample = mapDataLoss;
    case 'angle'
        mapDataTargetExample = mapDataAngle;
    case 'relaxance'
        mapDataTargetExample = mapDataRelaxance;
end
cid = find(strcmp({resultsStruct.(varNames{j}).clusterData.clusterVar}, clusterTarget));
mapDataCluster2DExample = resultsStruct.(varNames{j}).clusterData(cid).clusterMap2D;
mapDataGlobalCluster2DExample = resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap2D;

fprintf('complete!\n');

%% Figure 2 Data Selections
exampleCellType = {'HEMa','HFF','A375P','A375M1','A375M2'};     % Cell Type to use
exampleCellDish = {2, 1, 1, 2 ,1};                              % Dish Number of Type
exampleCellNumber = {2, 1, 1, 1, 1};                            % Cell Number in Dish
evalPt = [5e2, 1e3, 5e3, 10e3, 25e3, 50e3];                     % Frequency for evaluating maps
dataOrder = {'topo','storage','loss','angle'};                  % Order (L > R, Top > Bottom) to plot observables in quad

%% Figure 2
% For this figure, we create compact grids of viscoelastic observables for
% each cell type at a variety of frequencies. There is quite a lot of data,
% so runtimes for this section may be large.

fprintf('\nGenerating Figure 2...');
for i_file = 1:numel(exampleCellType)

    % Begin by finding the file we need for this figure and loading the data
    % Search recursively (using "**") for our zTransform results files
    dataFile = dir(fullfile(highResPath, '**',...
        ['*' num2str(exampleCellType{i_file}) 'Dish' num2str(exampleCellDish{i_file})...
        '*Cell' num2str(exampleCellNumber{i_file}) '_*Results*zTransform*.mat']));

    if isempty(dataFile)
        error('The directory you selected does not contain the Z-Transform QI map indicated for Figure 1. Please verify your MapResults file is in that directory and the filename contains "zTransform".');
    end

    resultsStruct = load([dataFile.folder filesep dataFile.name],'-mat');
    varNames = fields(resultsStruct);

    % If there are any other analysis results...skip them! This is
    % primarily for compatibility later, in case a new analysis type is
    % proposed in the future.    
    j = find(contains(varNames,'zTransform'));

    % Extract the size of the map from our results file
    if isfield(resultsStruct.(varNames{j}),'mapSize')
        mapSize = resultsStruct.(varNames{j}).mapSize;
    else
        mapSize = [128 128];
    end
    
    % Grab some relevant settings
    try
        correctTilt = resultsStruct.(varNames{j}).correctTilt;
        zeroSubstrate = resultsStruct.(varNames{j}).zeroSubstrate;
        optimizeFlattening = resultsStruct.(varNames{j}).optimizeFlattening;
    catch
        warning('The variables correctTilt, zeroSubstrate, and optimizeFlattening were not found in the results file. Using hard-coded settings in generatePaperFigures() (all set to true).');
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
    [minHeight,~] = min(pixelHeightArray(pixelHeightArray>0));
    substrateCutoff = minHeight + trimHeight;
    pixelsToRemove = false(size(pixelHeightArray));
    pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;
    
    % Remove the pixels we want to keep from the list
    pixelSkip = 1:numel(pixelHeightArray);
    pixelSkip(~pixelsToRemove) = [];
    
    % Make blank map data
    heightImg = zeros(flip(mapSize));
    heightImg = rot90(reshape(pixelHeightArray,mapSize),1);

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
        
        % Starting position for the map
        xc = 1;
        yc = 0;

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
                mapDataInd(idx_pixel) = interp1(t_t,h_t,1./freqList(k_freq),'makima',...
                    1e-12);
%                 mapDataHeight(idx_pixel) = pixelHeightArray(idx_pixel);
                mapDataHeight(idx_pixel) = heightImg(mapSize(2)-yc,xc);

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
            'TileSpacing', 'none');
%         pbaspect([1 1 1])

        for i_plot = 1:4

            ax = nexttile;
            plotTarget = dataOrder{i_plot};

            switch plotTarget
                case 'topo'
                    mapData = mapDataHeight;
                    plotTitle = 'Topography';
                    saveLabel = 'Topography';
                    
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
            grid on
            box(ax,boxSetting)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
%             title(plotTitle,'FontSize',largeFontSize)
%             xlabel(xlab,'FontSize',mediumFontSize)
%             ylabel(ylab,'FontSize',mediumFontSize)
            set(gca,'YTickLabel',[],'XTickLabel',[])
            xlim(xlims)
            ylim(ylims)
%             axis('equal')
            
            switch plotTarget
                case 'topo'
                    caxis(ax,[0 zMax]);
                case 'ind'
                    caxis(ax,[0 climInd]);
                case 'storage'
                    caxis(ax,[0 climMax]);
                case 'loss'
                    caxis(ax,[0 climMax]);
                case 'angle'
                    caxis(ax,[0 90]);
                case 'relaxance'
                    caxis(ax,[0 climMax]);
            end
            
            if showScaleBar && exist('XA','var') && exist('YA','var')
                barsize = 10^ceil(log10(max(xlims))-1);
                xbar = [(1/20)*xlims(2) (1/20)*xlims(2)+barsize];
                ybar = ones(size(xbar))*(1/10)*xlims(2);
                zbar = ones(size(ybar))*scaleBarMax;
                plot3(xbar,ybar,zbar,'-k','LineWidth',5)
                text(xbar(1),ybar(1),zbar(1),sprintf('%d \\mum',barsize),...
                    'HorizontalAlignment','left',...
                    'VerticalAlignment','top',...
                    'FontSize',tinyFontSize);

                xlim(xlims)
                ylim(xlims)
        %         zlim([0 scaleBarMax])
            end
            
            hold off
            view(2)
            
            if any(diff(yticks) ~= barsize)
                temp = yticks;
                yticks(ax,temp(1):barsize:temp(end))
            end
            if any(diff(xticks) ~= barsize)
                temp = xticks;
                xticks(ax,temp(1):barsize:temp(end))
            end
            
            pbaspect([1 1 1])
            drawnow
            
        end
        
        % Create our filename
        plotFile = [savePath filesep 'Fig2-QuadPlot-' exampleCellType{i_file} '-' num2str(freqList(k_freq)) 'Hz'];
        
        % Using exportgraphics() for Higher Quality
        saveas(mapPlotWindow,[plotFile '.fig'])
        exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
        exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);
    
    end
        
%     % Plot our topography and save it
%     % Clear old figures if they exist
%     if ~exist('mapPlotWindow','var')
%         mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
%     else
%         try
%             figure(mapPlotWindow)
%             clf
%             mapPlotWindow.Position = [figX figY figWid figHeight];
%         catch
%             mapPlotWindow = figure('Position',[figX figY figWid figHeight]);
%         end
%     end
% 
%     plotTitle = 'Topography';
%     saveLabel = 'Topography';
% 
%     surf(XPlot,YPlot,mapDataHeight,mapDataHeight,'EdgeColor','interp')
%     ax = gca;
%     colormap(ax,mapColorName)
%     hold on
%     grid on
%     box(ax,boxSetting)
%     set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
%     title(plotTitle,'FontSize',largeFontSize)
%     if showLabels
%         xlabel(xlab,'FontSize',mediumFontSize)
%         ylabel(ylab,'FontSize',mediumFontSize)
%     else
%         set(gca,'YTickLabel',[],'XTickLabel',[])
%     end
%     xlim(xlims)
%     ylim(ylims)
%     cb = colorbar;
%     caxis([0 zMax]);
%     temp = (cb.Ticks' ./ 1e-6);
%     for ii = 1:numel(temp)
%        cb.TickLabels{ii} = sprintf('%g \\mum',temp(ii));
%     end
%     hold off
% 
%     view(2)
%     pbaspect([1 1 1])
%     drawnow
% 
%     % Create our filename
%     plotFile = [savePath filesep 'Fig2-' saveLabel '-' exampleCellType{i_file}];
% 
%     % Using exportgraphics() for Higher Quality
%     saveas(mapPlotWindow,[plotFile '.fig'])
%     exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
%     exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);
    
    clearvars resultsStruct

end

fprintf('complete!\n');

%% Figure 3 Settings
colorList = turbo(10);        % Colors for Plotting

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
for i_plot = 1:4

    figure(mapPlotWindow)
    clf

    switch i_plot
        case 1
            mapData = mapDataHeightExample;
            plotTitle = 'Topography';
            saveLabel = 'Topography';
            
        case 2
            mapData = mapDataTargetExample;
            switch clusterTarget
                case 'force'
                    plotTitle = 'Force';
                    saveLabel = 'Force';
                case 'indentation'
                    plotTitle = 'Indentation';
                    saveLabel = 'Indentation';
                case 'storage'
                    plotTitle = 'Storage Modulus';
                    saveLabel = 'Storage';
                case 'loss'
                    plotTitle = 'Loss Modulus';
                    saveLabel = 'Loss';
                case 'angle'
                    plotTitle = 'Loss Angle';
                    saveLabel = 'Angle';
                case 'relaxance'
                    mapData = abs(mapDataTargetExample);
                    plotTitle = 'Relaxance';
                    saveLabel = 'Relaxance';
            end

        case 3
            % We want to order the clusters for this map according to
            % increasing magnitude (softest to largest) so all of our
            % clusters are aligned between the maps.
            uniqueBins = unique(mapDataCluster2DExample(~isnan(mapDataCluster2DExample)));
            medVal = NaN(size(uniqueBins));
            for k_align = 1:length(uniqueBins)
                % Grab our bin values at the frequency we will be plotting!
                % Then, we can find the median of that bin and use it to
                % sort our bin labels in increasing order.
                tempData = mapDataCluster2DExample(mapDataCluster2DExample == uniqueBins(k_align));
                medVal(k_align) = median(tempData,'omitnan');
            end
            [~,idold] = sort(medVal);
            idnew = uniqueBins(idold);
            clusterDataTemp = NaN(size(mapDataCluster2DExample));
            for k_align = 1:length(idold)
                % Reorganize our data according to the new order
                clusterDataTemp(mapDataCluster2DExample == idold(k_align)) = idnew(k_align);
            end
            mapData = clusterDataTemp;
%             mapData = mapDataCluster2DExample;
            plotTitle = 'Locally-Observed Clusters';
            saveLabel = 'SoloCluster';

        case 4
            % We want to order the clusters for this map according to
            % increasing magnitude (softest to largest) so all of our
            % clusters are aligned between the maps.
            uniqueBins = unique(mapDataGlobalCluster2DExample(~isnan(mapDataGlobalCluster2DExample)));
            medVal = NaN(size(uniqueBins));
            for k_align = 1:length(uniqueBins)
                % Grab our bin values at the frequency we will be plotting!
                % Then, we can find the median of that bin and use it to
                % sort our bin labels in increasing order.
                tempData = mapDataGlobalCluster2DExample(mapDataGlobalCluster2DExample == uniqueBins(k_align));
                medVal(k_align) = median(tempData,'omitnan');
            end
            [~,idold] = sort(medVal);
            idnew = uniqueBins(idold);
            clusterDataTemp = NaN(size(mapDataGlobalCluster2DExample));
            for k_align = 1:length(idold)
                % Reorganize our data according to the new order
                clusterDataTemp(mapDataGlobalCluster2DExample == idold(k_align)) = idnew(k_align);
            end
            mapData = clusterDataTemp;
%             mapData = mapDataGlobalCluster2DExample;
            plotTitle = 'Globally-Observed Clusters';
            saveLabel = 'GlobalCluster';

    end

    surf(XPlotExample,YPlotExample,mapDataHeightExample,...
        mapData,'EdgeColor','interp')
    ax = gca;
    if any(i_plot == [3,4])
        colormap(ax,colorList)
    else
        colormap(ax,mapColorName)
    end
    hold on
    grid on
    box(ax,boxSetting)
    set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
    cb = colorbar;
    set( findall(gcf, '-property', 'fontsize'), 'fontsize', smallFontSize)
    if showLabels
        title(plotTitle,'FontSize',largeFontSize)
        xlabel(xlabExample,'FontSize',mediumFontSize)
        ylabel(ylabExample,'FontSize',mediumFontSize)
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end
    xlim(xlimsExample)
    ylim(ylimsExample)
    switch i_plot
        case 1
            caxis([0 zMax]);
            temp = (cb.Ticks' ./ 1e-6);
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%g',temp(ii)); % \\mum
            end
            scaleBarMax = zMax;
            
        case 2
            switch clusterTarget
                case 'force'
                    caxis([0 climForce]);
                    temp = (cb.Ticks' ./ 1e-12);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%g',temp(ii)); % pN
                    end
                    scaleBarMax = climForce;
                case 'indentation'
                    caxis([0 climInd]);
                    temp = (cb.Ticks' ./ 1e-9);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%g',temp(ii)); % nm
                    end
                    scaleBarMax = climInd;
                case 'storage'
                    caxis([0 climMax]);
                    temp = (cb.Ticks' .* 1e-3);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%d',round(temp(ii))); % kPa
                    end
                    scaleBarMax = climMax;
                case 'loss'
                    caxis([0 climMax]);
                    temp = (cb.Ticks' .* 1e-3);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%d',round(temp(ii))); % kPa
                    end
                    scaleBarMax = climMax;
                case 'angle'
                    cb.Ruler.TickLabelFormat='%g'; % Deg
                    caxis([0 90]);
                    cb.Ticks = [0 30 60 90];
                    scaleBarMax = 90;
                case 'relaxance'
                    caxis([0 climMax]);
                    temp = (cb.Ticks' .* 1e-3);
                    for ii = 1:numel(temp)
                       cb.TickLabels{ii} = sprintf('%d',round(temp(ii))); % kPa
                    end
                    scaleBarMax = climMax;
            end
                    
        case {3,4}
%             caxis([1 max(unique(mapData(~isnan(mapData))))]);
            caxis([0 size(colorList,1)]);
            cb.Ticks = 0:size(colorList,1);
            temp = cb.Ticks';
            for ii = 1:numel(temp)
               cb.TickLabels{ii} = sprintf('%d',round(temp(ii)));
            end
            
            % Remove the 0-point
            cb.TickLabels{1} = [];
            
            % Shift the labels
            cb.Ticks = cb.Ticks - 0.5;

    end
        
    view(2)
    
    if showScaleBar && exist('XAExample','var') && exist('YAExample','var')
        
        xlims = [0 max(scanSizeExample)./1e-6];
        
        barsize = 10^ceil(log10(max(xlims))-1);
        xbar = [(1/20)*xlims(2) (1/20)*xlims(2)+barsize];
        ybar = ones(size(xbar))*(1/10)*xlims(2);
        zbar = ones(size(ybar))*scaleBarMax;
        plot3(xbar,ybar,zbar,'-k','LineWidth',5)
        text(xbar(1),ybar(1),zbar(1),sprintf('%d \\mum',barsize),...
            'HorizontalAlignment','left',...
            'VerticalAlignment','top',...
            'FontSize',smallFontSize);

        xlim(xlims)
        ylim(xlims)
%         zlim([0 scaleBarMax])
    end

    hold off
    view(2)

    if any(diff(yticks) ~= barsize)
        temp = yticks;
        yticks(ax,temp(1):barsize:temp(end))
    end
    if any(diff(xticks) ~= barsize)
        temp = xticks;
        xticks(ax,temp(1):barsize:temp(end))
    end
    
    pbaspect([1 1 1])
    drawnow

    % Create our filename
    plotFile = [savePath filesep 'Fig3-' saveLabel];
    
    % Using exportgraphics() for Higher Quality
    saveas(mapPlotWindow,[plotFile '.fig'])
    exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
    exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

end

fprintf('complete!\n');

%% Figure 4 Data Selections
clusterCellTypes = {'HEMa','HFF','A375P','A375M1','A375M2'};    % Cell Type to use
evalPt = [1e3 10e3];                                                   % Frequency for evaluating maps
% nBins = 5000;                                                   % Number of bins for the histogram
climMaxHist = climMax;                                          % Pa, maximum stress to show on main histogram
climMaxHistAll = {[0 20e3],[90e3 110e3],[180e3 200e3]};         % Pa, limits of subplots to generate for each histogram            

% % Grab bin designations
% [~,sheet_names] = xlsfinfo('binLabels.xlsx');
% binLabels = cell(numel(sheet_names),1);
% for k = 1:numel(sheet_names)
%   [~,~,binLabels{k}] = readtable('binLabels.xlsx',sheet_names{k});
% end

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
    Files = dir(fullfile(clusterPath, '**',...
        ['*' clusterCellTypes{i_dir} '*Results*zTransform*.mat']));

    if isempty(Files)
        error('The directory you selected and its subdirectories do not contain a Z-Transform QI map for cell type No. %d. Please verify your results file is in or under that directory and the filename contains "zTransform".',i_dir);
    end

    % Since we are analyzing many maps, we must create our frequency
    % array ahead of time.
    minFreq = Inf;
    maxFreq = 0;
    mapsizes = NaN(numel(Files),1);
    
    % Prepare our index records
    ai = NaN(numel(Files),1);
    bi = NaN(numel(Files),1);
    for j_dir = 1:length(Files)
        resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
        mapsizes(j_dir) = numel(resultsStruct.zTransformAnalysis.frequencyMap);
        clearvars resultsStruct
    end
    
    for j_dir = 1:numel(mapsizes)
        if j_dir == 1
            ai(j_dir) = 1;
        else
            ai(j_dir) = bi(j_dir-1) + 1;
        end
        bi(j_dir) = ai(j_dir) + mapsizes(j_dir) - 1;
    end

    for j_dir = 1:length(Files)
        
        vars = whos('-file',[Files(j_dir).folder filesep Files(j_dir).name]);
        varNames = {vars.name};
        
        % If there are any other analysis results...skip them! This is
        % primarily for compatibility later, in case a new analysis type is
        % proposed in the future.    
        j = find(contains(varNames,'zTransform'));
        
        % Load the mat file
        resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
        freqMapTemp = resultsStruct.(varNames{j}).frequencyMap;
        timesCellTemp = resultsStruct.(varNames{j}).ViscoClass.times_cell;
        clearvars resultsStruct

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
    
%     clusteringDataGlobal = NaN(bi(end),numel(magList));
    XDataGlobal = NaN(bi(end),length(evalPt));
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
        j = find(contains(varNames,'zTransform'));

        if isfield(resultsStruct.(varNames{j}),'mapSize')
            mapSize = resultsStruct.(varNames{j}).mapSize;
        else
            mapSize = [128 128];
        end

        % Grab some relevant settings
        try
            correctTilt = resultsStruct.(varNames{j}).correctTilt;
            hideSubstrate = resultsStruct.(varNames{j}).hideSubstrate;
            zeroSubstrate = resultsStruct.(varNames{j}).zeroSubstrate;
            optimizeFlattening = resultsStruct.(varNames{j}).optimizeFlattening;
        catch
            hideSubstrate = true;
            warning('The variables hideSubstrate, correctTilt, zeroSubstrate, and optimizeFlattening were not found in the results file. Using hard-coded settings in generatePaperFigures() (all set to true).');
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

        pixelLog = NaN(numel(pixelHeightArray),2);
        histXData = NaN(numel(pixelHeightArray),length(evalPt));

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

                for i_eval = 1:length(evalPt)
                    obsOut = interp1(xinterp,clusterInterp,evalPt(i_eval),'makima',...
                        NaN);
                    histXData(k_pixels,i_eval) = obsOut;
                end
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
        histYData = reshape(resultsStruct.(varNames{j}).clusterData(cid).globalClusterMap2D,bi(j_dir)-(ai(j_dir)-1),1);

        
%         for i_eval = 1:length(evalPt)
            
            % We want to order the clusters for this map according to
            % increasing magnitude (softest to largest) so all of our
            % clusters are aligned between the maps.
            uniqueBins = unique(histYData(~isnan(histYData)));
            medVal = NaN(size(uniqueBins));
            for k_align = 1:length(uniqueBins)
                % Grab our bin values at the frequency we will be plotting!
                % Then, we can find the median of that bin and use it to
                % sort our bin labels in increasing order. For more than 1
                % evalPt, we baseline on the first frequency.
                tempData = histXData(histYData == uniqueBins(k_align),1);
%                 tempData = histXData(histYData == uniqueBins(k_align),i_eval);
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
            
%         end

        clearvars resultsStruct mapSize histXData histYData histYDataTemp ...
            pixelLog

    end
    fprintf('complete! The %s dataset contains %d pixels.\n',...
        clusterCellTypes{i_dir},size(XDataGlobal,1));

    % Now, go through and save the results to each output file and make
    % our plots!
    fprintf('\nGenerating our histogram plots...');

    % Loop through all the eval points
    for i_eval = 1:length(evalPt)
    
        % Clear old figures if they exist
        if ~exist('mapPlotWindow','var')
            mapPlotWindow = figure('Position',[figX figY figWid*2 figHeight*0.6]);
        else
            try
                figure(mapPlotWindow)
                clf
                mapPlotWindow.Position = [figX figY figWid*2 figHeight*0.6];
            catch
                mapPlotWindow = figure('Position',[figX figY figWid*2 figHeight*0.6]);
            end
        end

        % Now, begin making our histogram
        uniqueBins = unique(YDataGlobal(~isnan(YDataGlobal)));
        hi = gobjects(numel(uniqueBins),1);
        ax = gca;

        % Determine plot order
        pixCount = zeros(size(uniqueBins));
        for i_plot = 1:numel(uniqueBins)
            % We don't want to hide any smaller distributions! So we have to
            % plot from largest-to-smallest histogram (like stacking blocks).
            % This way, we hope, none of the distributions are hidden but they
            % are still labeled correctly.
            pixCount(i_plot) = sum(YDataGlobal == uniqueBins(i_plot));
        end
        [~,sortid] = sort(pixCount,'descend');
        uniqueBins = uniqueBins(sortid);

        % Create and trim our data
        histDatasetTemp = XDataGlobal(:,i_eval);
        switch clusterTarget
            case 'force'
                histMax = climForce;
            case 'indentation'
                histMax = climInd;
            case {'storage','loss','relaxance'}
                histMax = climMaxHist;
            case 'angle'
                histMax = 90;
        end
        histDatasetTemp(histDatasetTemp > histMax | histDatasetTemp < 0) = [];
        [~,binEdges] = histcounts(histDatasetTemp,...
            'BinMethod','scott');
%         [~,binEdges] = histcounts(histDatasetTemp,20); % Manual number of bins
        clearvars histDatasetTemp

        hold on
        grid on
        for i_plot = 1:numel(uniqueBins)

            if exist('binLabels','var')
                % Find out which color to use
                % THIS CODE NEEDS TO BE TESTED!!!

                sheetidx = find(strcmpi(sheet_names,clusterCellTypes{i_dir})); % Pick the right sheet for this type of cell
                temp = find(binLabels{1}{:,2}==binLabels{sheetidx}{uniqueBins(i_plot),2}); % Find the name of this bin in the first sheet (roster)
                colorID = binLabels{1}{temp,1}; % Pick the "bin number" associated with that name
            else
                colorID = uniqueBins(i_plot);
            end

            % Create and trim our data
            histDataset = XDataGlobal(YDataGlobal == uniqueBins(i_plot),i_eval);
            histDataset(histDataset > histMax | histDataset < 0) = [];
    %         dBin = (max(histDataset) - min(histDataset))/numel(binCounts);
    %         nBins = numel(binCounts); % Adjust the number using the resolution predicted by histcounts
    %         binEdges = linspace(0,climMaxHist,nBins+1);
    %         binEdges = linspace(0,climMax,nBins+1);

            % Create our histogram
            hi(i_plot) = histogram(histDataset,...
                'BinEdges', binEdges,...
                'FaceColor', colorList(colorID,:),...
                'FaceAlpha', 0.9,...
                'EdgeColor', colorList(colorID,:),...
                'EdgeAlpha', 0.5);
    %             'FaceColor', colorList(colorID,:),...
    %             'FaceAlpha', 0.5,...
    %             'EdgeColor', [0 0 0],...
    %             'EdgeAlpha', 0.5);

    %         % Plot a normal distribution over the top
    %         [mu,sigma] = normfit(hi(i_plot).Values);
    %         xvals = binEdges(1:end-1) + diff(binEdges)/2;
    %         yvals = max(hi(i_plot).Values,[],'all').*exp(-(hi(i_plot).Values-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
    %         plot(xvals,yvals,'Color',colorList(i_plot,:),'LineWidth',1,'HandleVisibility','off')
        end
        box(ax,boxSetting)
        set(gca,'TickLength',[0.01 0.01],'LineWidth',2)
        set( findall(gcf, '-property', 'fontsize'), 'fontsize', smallFontSize)

        switch clusterTarget
            case 'force'
                xlim([0 climForce])
                xlab = sprintf('Force [pN]');
                temp = (xticks' ./ 1e-12);
                TickLabels = cell(size(temp));
                for ii = 1:numel(temp)
                   TickLabels{ii} = sprintf('%d',round(temp(ii)));
                end
                xticklabels(TickLabels)
            case 'indentation'
                xlim([0 climInd])
                xlab = sprintf('Indentation [\\mum]');
                temp = (xticks' ./ 1e-6);
                TickLabels = cell(size(temp));
                for ii = 1:numel(temp)
                   TickLabels{ii} = sprintf('%g',temp(ii));
                end
                xticklabels(TickLabels)
            case 'storage'
                xlim([0 climMaxHist])
                xlab = sprintf('Storage Modulus [kPa]');
                temp = (xticks' .* 1e-3);
                TickLabels = cell(size(temp));
                for ii = 1:numel(temp)
                   TickLabels{ii} = sprintf('%d',round(temp(ii)));
                end
                xticklabels(TickLabels)
            case 'loss'
                xlim([0 climMaxHist])
                xlab = sprintf('Loss Modulus [kPa]');
                temp = (xticks' .* 1e-3);
                TickLabels = cell(size(temp));
                for ii = 1:numel(temp)
                   TickLabels{ii} = sprintf('%d',round(temp(ii)));
                end
                xticklabels(TickLabels)
            case 'angle'
                xlim([0 90])
                xlab = sprintf('Loss Angle [deg]');
            case 'relaxance'
                xlim([0 climMaxHist])
                xlab = sprintf('Relaxance [kPa]');
                temp = (xticks' .* 1e-3);
                TickLabels = cell(size(temp));
                for ii = 1:numel(temp)
                   TickLabels{ii} = sprintf('%d',round(temp(ii)));
                end
                xticklabels(TickLabels)
        end
%         ylabel('Number of Pixels', 'FontSize', mediumFontSize)

        if i_dir == length(clusterCellTypes)
            xlabel(xlab, 'FontSize', mediumFontSize)
        else
            set(gca,'Xticklabel',[])
        end

        legendEntries = sprintfc('%d',uniqueBins);
        [~,sortidx] = sort(uniqueBins,'ascend');
        legend(hi(sortidx),legendEntries(sortidx),'location','eastoutside')

        hold off

        % Create our filename
        plotFile = [savePath filesep 'Fig4-' clusterCellTypes{i_dir} sprintf('-%dHz',evalPt(i_eval))];

        % Using exportgraphics() for Higher Quality
        saveas(mapPlotWindow,[plotFile '.fig'])
        exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
        exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

        % Now, do the same thing for all the subplots
        for i_subs = 1:length(climMaxHistAll)

            % Clear old figures if they exist
            if ~exist('mapPlotWindow','var')
                mapPlotWindow = figure('Position',[figX figY figWid*1.4 figHeight*0.8]);
            else
                try
                    figure(mapPlotWindow)
                    clf
                    mapPlotWindow.Position = [figX figY figWid*1.4 figHeight*0.8];
                catch
                    mapPlotWindow = figure('Position',[figX figY figWid*1.4 figHeight*0.8]);
                end
            end

            % Now, begin making our histogram
            uniqueBins = unique(YDataGlobal(~isnan(YDataGlobal)));
            hi = gobjects(numel(uniqueBins),1);
            ax = gca;

            % Create and trim our data
            histDatasetTemp = XDataGlobal(:,i_eval);
            histBinTemp = YDataGlobal;
            histMax = climMaxHistAll{i_subs};
            histBinTemp(histDatasetTemp < histMax(1) | histDatasetTemp > histMax(2)) = [];
            histDatasetTemp(histDatasetTemp < histMax(1) | histDatasetTemp > histMax(2)) = [];
            [Nbin,binEdges] = histcounts(histDatasetTemp,...
                'BinMethod','scott');
            if numel(Nbin) < 20
                [~,binEdges] = histcounts(histDatasetTemp,20); % Manual number of bins
            end
            
            % Determine plot order
            pixCount = zeros(size(uniqueBins));
            for i_plot = 1:numel(uniqueBins)
                % We don't want to hide any smaller distributions! So we have to
                % plot from largest-to-smallest histogram (like stacking blocks).
                % This way, we hope, none of the distributions are hidden but they
                % are still labeled correctly.
%                 pixCount(i_plot) = sum(YDataGlobal == uniqueBins(i_plot));
                pixCount(i_plot) = sum(histBinTemp == uniqueBins(i_plot));
            end
            [~,sortid] = sort(pixCount,'descend');
            uniqueBins = uniqueBins(sortid);
            
            clearvars histDatasetTemp

            hold on
            grid on
            for i_plot = 1:numel(uniqueBins)

                % Create and trim our data
                histDataset = XDataGlobal(YDataGlobal == uniqueBins(i_plot),i_eval);
                histDataset(histDataset < histMax(1) | histDataset > histMax(2)) = [];

                % Create our histogram
                hi(i_plot) = histogram(histDataset,...
                    'BinEdges', binEdges,...
                    'FaceColor', colorList(uniqueBins(i_plot),:),...
                    'FaceAlpha', 0.9,...
                    'EdgeColor', colorList(uniqueBins(i_plot),:),...
                    'EdgeAlpha', 0.5);

            end
            box(ax,boxSetting)
            set(gca,'TickLength',[0.01 0.01],'LineWidth',2)
            set( findall(gcf, '-property', 'fontsize'), 'fontsize', mediumFontSize)

            switch clusterTarget
                case 'force'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Force [pN]');
                    temp = (xticks' ./ 1e-12);
                    TickLabels = cell(size(temp));
                    for ii = 1:numel(temp)
                       TickLabels{ii} = sprintf('%d',round(temp(ii)));
                    end
                    xticklabels(TickLabels)
                case 'indentation'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Indentation [\\mum]');
                    temp = (xticks' ./ 1e-6);
                    TickLabels = cell(size(temp));
                    for ii = 1:numel(temp)
                       TickLabels{ii} = sprintf('%g',temp(ii));
                    end
                    xticklabels(TickLabels)
                case 'storage'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Storage Modulus [kPa]');
                    temp = (xticks' .* 1e-3);
                    TickLabels = cell(size(temp));
                    for ii = 1:numel(temp)
                       TickLabels{ii} = sprintf('%d',round(temp(ii)));
                    end
                    xticklabels(TickLabels)
                case 'loss'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Loss Modulus [kPa]');
                    temp = (xticks' .* 1e-3);
                    TickLabels = cell(size(temp));
                    for ii = 1:numel(temp)
                       TickLabels{ii} = sprintf('%d',round(temp(ii)));
                    end
                    xticklabels(TickLabels)
                case 'angle'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Loss Angle [deg]');
                case 'relaxance'
                    xlim([histMax(1) histMax(2)])
                    xlab = sprintf('Relaxance [kPa]');
                    temp = (xticks' .* 1e-3);
                    TickLabels = cell(size(temp));
                    for ii = 1:numel(temp)
                       TickLabels{ii} = sprintf('%d',round(temp(ii)));
                    end
                    xticklabels(TickLabels)
            end
%             ylabel('Number of Pixels', 'FontSize', largeFontSize)
            xlabel(xlab, 'FontSize', largeFontSize)

            % Remove 0-mark yticklabel
            ax = gca;
            ytickstr = get(ax, 'YTickLabel');
            idx = strcmpi(ytickstr,'0');
            ytickstr{idx} = '';   % needs to exist but make it empty
            set(ax, 'YTickLabel', ytickstr);
            clearvars ax
            
            hold off

            % Create our filename
            plotFile = [savePath filesep 'Fig4-' clusterCellTypes{i_dir} sprintf('-%dHz',evalPt(i_eval)) '-SubplotRange' sprintf('%d',i_subs)];

            % Using exportgraphics() for Higher Quality
            saveas(mapPlotWindow,[plotFile '.fig'])
            exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
            exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

        end
        
    end
    
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