function [] = singleCurveClustering(originalPath,N_workers,clusterTarget,varargin)
%SINGLECURVECLUSTERING Perform K-Medoids Clustering of a group of SFS
%Experiments
%   This function takes in a path argument, the target variable for
%   clustering with k-medoids, and a variable number of arguments.

%   Note that evalPt is the frequency (in Hz) to use in our analysis
%   for plotting the frequency-dependent properties from the data. The
%   full time series are available for each curve and are used in the
%   clustering process, and evalPt is exclusively for visualization
%   purposes.

% User-Defined Settings
logSteps = true;
n_reps = 10;                % number of clustering replicates
maxK = 10;                  % Max number of cluster bins
if nargin > 1
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        logSteps = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        n_reps = varargin{i};                        
                    end
                case 3
                    if ~isempty(varargin{i})
                        maxK = varargin{i};                        
                    end
                otherwise
                    fprintf('Passed additional parameters to singleCurveClustering() which were not used.');
            end
        end
    end
end

% Permanent Settings
figX = 0;
figY = 0;
maxCol = 3;
maxwid = get(0,'screensize');
maxwid = maxwid(3);
titleFontSize = 14;
axisFontSize = 14;
n_steps = 10; % frames, number of frames per order of magnitude

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
        
        % Load all of the force curves available in our directory!
        % Settings for loading the data
        loadDataSettings  = struct();

        % Required Settings:
        loadDataSettings.includeRetract = false;         % Don't include data from the retract curve
        loadDataSettings.filterType = 'none';            % Choose the filter used to smooth data
        loadDataSettings.findRep = 'legendre';           % Search direction for the repulsive region
        loadDataSettings.removeNegatives = true;         % Remove negative values in the force
        loadDataSettings.createAverage = false;          % Create averaged rows at the END of the datastruct

        % Conditional Settings (depending on filter):
        loadDataSettings.N = 2;                          % Order of Butterworth filter (if used)
        loadDataSettings.cutoff_Hz = 5000;               % Cutoff frequency of Butterworth (if used)

        % Load the AFM Data
        dataStruct = LoadAFMData(path,loadDataSettings);
        
        if isempty(dataStruct)
            error('The directory you selected does not contain any usable data. Please verify your files are in that directory and the filetypes are valid.');
        end

        fileLabels = cell(numel(dataStruct),1);
        for j = 1:numel(dataStruct)
            temp = strsplit(path,{filesep});
            fileLabels{j} = strjoin({temp{end} num2str(j)},'-');
        end
            
        minTime = 0;
        maxTime = Inf;
        dtlist = NaN(numel(dataStruct),1);
        for i = 1:numel(dataStruct)
            dtlist(i) = dataStruct(i).dt;
            tempt = dataStruct(i).t_r;
            tempt = tempt(tempt>0);
            [temp,~] = min(tempt,[],'omitnan');
            if temp > minTime
                minTime = temp;
            end
            tempt = dataStruct(i).t_r;
            tempt = tempt(tempt>0);
            [temp,~] = max(tempt,[],'omitnan');
            if temp < maxTime
                maxTime = temp;
            end
        end
        
        % Count orders of 10
        temp = minTime;
        tempt = minTime;
        tempmax = 10^ceil(log10(minTime));
        magList = [];
        if ~logSteps
            dt = min(dtlist,[],'all'); % Hz, step size between frames, if discrete
            while temp < maxTime
                tempt = 10.^( ( log10(temp) ) );
                magList = horzcat(magList,tempt);
                temp = temp + dt;
            end
        else
            while temp < maxTime
                if temp == minTime
                    dt = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_steps);
                end

                if temp >= tempmax
                    tempmax = temp*10;
                    dt = ((10^(ceil(log10(temp)))-10^(floor(log10(temp))))/n_steps);
                end

                tempt = 10.^( ( log10(temp) ) );
                magList = horzcat(magList,tempt);
                temp = temp + dt;
            end
        end
        magList = unique(magList);
        
        switch clusterTarget
            case 'force'
                saveLabel = 'Force';
            case 'indentation'
                saveLabel = 'Ind';
        end

        plotFile = [path filesep 'SingleCurvesClustering-BinPlots-' saveLabel];
        clusteringData = NaN(numel(dataStruct),numel(magList));
        
        parfor k_curve = 1:numel(dataStruct)

            % Resample to known array of times
            t_t = dataStruct(k_curve).t_r;
            ids = ((t_t >= min(magList)) & (t_t <= max(magList)));

            if sum(ids,'all') < 2
                % Skip datasets with only one datapoint
                continue;
            end

            h_t = dataStruct(k_curve).h_r;
            F_t = dataStruct(k_curve).F_r;
                    
            clusterInterp = [];
            switch clusterTarget
                case 'force'
                    clusterInterp = F_t;
                case 'indentation'
                    clusterInterp = h_t;
            end

            try
                % Resample to known array of frequencies
                obsOut = interp1(t_t,clusterInterp,magList,'makima',...
                    NaN);
                clusteringData(k_curve,:) = obsOut;
            catch
                % Do nothing
            end
                    
        end
        
        opts = statset('UseParallel',parallelSet,...
            'MaxIter',1e2,...
            'Display','off');
        tempfunc = @(x,k) kmedoidsnan(x,k,'Options',opts,...
            'Distance',@dtwf,...
            'Replicates',n_reps);
        ids = all(isnan(clusteringData),2); % Find excluded pixels
        clusteringData(ids,:) = [];
        
        % Try all of the cluster configurations. We start with 3
        % clusters to account for any substrate pixels that were
        % not trimmed. This is not ideal, but unavoidable with
        % automatic large-scale analysis.
        eva = evalclusters(clusteringData,tempfunc,'CalinskiHarabasz',...
            'klist',2:min([maxK size(clusteringData,1)],[],'all'));

        fprintf('\nOptimal Number of Bins: %d\n\n',eva.OptimalK);

        idxK = tempfunc(clusteringData,eva.OptimalK);
        
        n_plots = eva.OptimalK;
        n_rows = ceil(n_plots/maxCol);
        n_cols = min([n_plots maxCol]);              
        mult = min([500 maxwid/n_cols]);
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
        
        tiledlayout(n_rows,n_cols, 'padding', 'none', ...
            'TileSpacing', 'tight', ...
            'OuterPosition', [0 0.15 1 0.85])

        % Loop through subplots
        for i_plot = 1:n_plots
        
            curveIDs = find(idxK==i_plot);
            n_curves = numel(curveIDs);
            cspecList = hsv(n_curves);
            
            ax = nexttile;

            hold on
            grid on
            box on
            title(sprintf('Cluster Bin %d',i_plot),'FontSize',titleFontSize)
            switch clusterTarget
                case 'force'
                    ylabel('Force [$$nN$$]','interpreter','latex','FontSize',axisFontSize)
                case 'indentation'
                    ylabel('Indentation [$$\mu m$$]','interpreter','latex','FontSize',axisFontSize)
            end
            xlabel('Time [$$s$$]','interpreter','latex','FontSize',axisFontSize)
            xlim([0 max(magList)])
            ylim([0 max(clusteringData,[],'all','omitnan')])
            
            % Loop through individual curves
            % Here, we use the data we actually used in clustering, NOT the
            % raw datasets. This gives us a direct view of the information
            % used to relate curves to one another, NOT what they're based
            % off of. If the curves are low quality, we need to change the
            % way we sample the curves above when generating the
            % clusteringData variable!
            for j_plot = 1:n_curves
                id = curveIDs(j_plot);
                xdata = magList;
                ydata = clusteringData(id,:);
                xdata(isnan(ydata)) = [];
                ydata(isnan(ydata)) = [];
                plot(ax,xdata,ydata,'LineStyle','-','Color',cspecList(j_plot,:))
            end
            
            % Updating the labels on our axes (aesthetic)
%             TickLabels = {};
%             switch clusterTarget
%                 case 'force'
%                     temp = (yticks' ./ 1e-6);
%                     for ii = 1:numel(temp)
%                        TickLabels{ii} = sprintf('%.2f \\muN',temp(ii));
%                     end
%                     yticklabels(ax,TickLabels)
%                 case 'indentation'
%                     temp = (yticks' ./ 1e-9);
%                     for ii = 1:numel(temp)
%                        TickLabels{ii} = sprintf('%.2f nm',temp(ii));
%                     end
%                     yticklabels(ax,TickLabels)
%             end
            
            % Force the plot to be square using pbaspect()
            pbaspect([1 1 1])
            
            hold off

        end
        
        saveas(mapPlotWindow,[plotFile '.fig'])
        print(mapPlotWindow,[plotFile '.png'],'-dpng','-r300');
        
        clusterResults = struct;
        
        % Begin creating our output file
        clusterResults.rawData = dataStruct;
        clusterResults.resampleTime = magList;
        clusterResults.clusteringData = clusteringData;
        clusterResults.fileLabels = fileLabels;
        
        % Create a placeholder struct where we store all of the
        % cluster results for each type of analysis. This is
        % critical, since it means we can run studies on
        % different clustering methods and compare results
        % later.
        temp = struct;
        for jj = 1:2
            switch jj
                case 1
                    temp(jj).clusterVar = 'force';

                case 2
                    temp(jj).clusterVar = 'indentation';                             
            end
            temp(jj).clusterBins = {};
            temp(jj).lastUpdate = '';
        end

        clusterResults.clusterData = temp;
        cid = find(strcmp({clusterResults.clusterData.clusterVar}, clusterTarget));
        clusterResults.clusterData(cid).clusterBins = num2cell(idxK');
        clusterResults.clusterData(cid).lastUpdate = datestr(now);
                
        % Save a results file
        temp = strsplit(path,{filesep});
        save([path filesep 'ClusteringResults-' temp{end} '-' saveLabel],'-struct','clusterResults','-v7.3');

    catch ERROR
        
        fprintf('ERROR Clustering Curves in Directory #%d of %d\n',i_dir,length(Folders));
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