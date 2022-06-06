clear all
close all
clc

% For a 128x map
pixelIDs = [500 1000 5000 8000 12500 14000 15025 15950];
savePath = uigetdir(pwd,'Select Output Dir');
evalPt = 1000;

%% Settings
hideSubstrate = true;
fillPixels = true;
logSteps = true;
showLabels = false;
showScaleBar = true;
clusterTarget = 'storage';
climMax = 2e5;

figX = 0;
figY = 0;
maxwid = get(0,'screensize');
maxheight = maxwid(4);
maxwid = maxwid(3);
figWid = maxwid;
figHeight = maxheight;
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
smallFontSize = 16;
mediumFontSize = 24;
largeFontSize = 32;

% Placeholders if our file doesn't have these settings
correctTilt = true;
zeroSubstrate = true;
optimizeFlattening = true;

%% Analyze
try
resultsStruct = load('D:\AFM Data\! HPC\ViscoFitZ-HPC\HighResViscoelasticityManuscript\4 - Second Final Data (Post-Global)\HFF\27-January-2022\Dish1\Cell4\HFFDish1_Cell4_MapResults_zTransform.mat','-mat');
catch
    [filename, pathname] = uigetfile('*.mat','Select a "MapResults" file');
    load([pathname filesep filename],'-mat');
end
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

        if ~any(ismember(k_pixels,pixelIDs))
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
mapDataHeight = rot90(reshape(pixelHeightArray,mapSize),1);
mapDataStorage = NaN(flip(mapSize));
mapDataLoss = NaN(flip(mapSize));
mapDataAngle = NaN(flip(mapSize));
mapDataRelaxance = NaN(flip(mapSize));

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
    
    mapDataHeight(idx_pixel) = heightImg(mapSize(2)-yc,xc);
    
    % Skip the pixel if we have decided to hide the substrate AND
    % this pixel is in our ignore list based on it's height
    if any(ismember(k_pixels,pixelIDs))
    
        % Make our observable plot
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
        
        tiledlayout(1,2)
        nexttile
        
        % height plot
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

        ax = gca;

        % We are hiding the colobars here.
        surf(XPlot,YPlot,mapDataHeight,mapDataHeight,'EdgeColor','interp')
        colormap(ax,mapColorName)
        hold on
        grid on
        box(ax,boxSetting)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
        set(gca,'YTickLabel',[],'XTickLabel',[])
        xlim(xlims)
        ylim(ylims)
        caxis(ax,[0 zMax]);

        if showScaleBar && exist('XA','var') && exist('YA','var')
            barsize = 10^ceil(log10(max(xlims))-1);
            xbar = [(1/20)*xlims(2) (1/20)*xlims(2)+barsize];
            ybar = ones(size(xbar))*(1/10)*xlims(2);
            zbar = ones(size(ybar))*zMax;
            plot3(xbar,ybar,zbar,'-w','LineWidth',5)
            text(xbar(1),ybar(1),zbar(1),sprintf('%d \\mum',barsize),...
                'HorizontalAlignment','left',...
                'VerticalAlignment','top',...
                'FontSize',tinyFontSize,...
                'Color','white');

            xlim(xlims)
            ylim(xlims)
        end
        
        % Highlight the pixel
        xpos = pixelLog(k_pixels,2);
        ypos = pixelLog(k_pixels,1);
        scatter3(XPlot(ypos,xpos),YPlot(ypos,xpos),mapDataHeight(ypos,xpos)*1.1,50,'rs','filled')
        
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

        % Plot Z-Transform viscoelastic observables example datasets
        % First, find ylim max because both axes have the same units
        nexttile
        
        Estorage_hz = modelStorage;
        Eloss_hz = modelLoss;
        tempMax = [1.1*max(Estorage_hz) 1.1*max(Eloss_hz)];
        pbaspect([1 1 1])

        % Now, make the plot
        yyaxis left
        ax = gca;
        plot(freq,Estorage_hz,'LineWidth',8)
        % plot(freq,smooth(Estorage_hz,floor(numel(Estorage_hz)*0.1)),'LineWidth',8)
        hold on
        grid on
        box(ax,boxSetting)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',2)
        set( findall(gca, '-property', 'fontsize'), 'fontsize', smallFontSize)
        xlabel('Frequency [$$kHz$$]','Interpreter','latex','FontSize',mediumFontSize)
        xlim([0 max(freq)])
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
        plot(freq,Eloss_hz,'LineWidth',8)
        % plot(freq,smooth(Eloss_hz,floor(numel(Eloss_hz)*0.1)),'LineWidth',5)
        hold on
        grid on
        box(ax,boxSetting)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',2)

        xlim([0 max(freq)])
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
        plotFile = [savePath filesep 'pixelobs-ZViscoModuli-Pixel' sprintf('%d',k_pixels)];
        saveas(mapPlotWindow,[plotFile '.fig'])
        exportgraphics(mapPlotWindow,[plotFile '.jpg'],'Resolution',300);
        exportgraphics(mapPlotWindow,[plotFile '.png'],'Resolution',300);

    end
    
    xc = xc + 1;

end
