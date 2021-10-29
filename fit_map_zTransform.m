function [] = fit_map_zTransform(originalPath,savePrependIn,N_workers,varargin)

if ischar(N_workers), N_workers = str2num(N_workers); end

correctTilt = false;
hideSubstrate = false;
if nargin > 3
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
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

format long
clc

% ======================== %
% Viscoelastic Parameter
% Extraction Script, Viscofit
% Z-Transform Version
%
% Created by cparvini
% ======================== %
% This file will perform the 
% required data analysis on 
% AFM SFS curves to extract 
% meaningful information about 
% how stiff the samples are on 
% specific timescales. It is an
% adaptation of the original code
% used specifically for compiling,
% as is necessary to use for 
% certain HPC facilities.
% ======================== %
% To run this, you need to send
% "savePrependIn", "originalPath",
% and "N_workers" from the SLURM 
% script. You can also pass two
% additional parameters that 
% control whether the substrate is
% excluded from analysis, and
% whether to fix any tilt issues
% with the QI map.

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

for i_dir = 1:length(Folders)
    
    % Start with a clean slate
    clearvars -except i_dir Folders originalPath savePrependIn N_workers correctTilt hideSubstrate
    close all
    clc
    fprintf('Analyzing Directory #%d of %d\n',i_dir,length(Folders));
    
    % Give the files a unique prepend, if necessary
    if length(Folders) > 1
        temp = strsplit(Folders{i_dir},filesep);
        savePrepend = [savePrependIn '_' temp{end}];
    else
        savePrepend = savePrependIn;
    end
    
    % Clean up the workers (memory management)
    if ~isempty(gcp('nocreate'))
        poolobj = gcp('nocreate');
        parfevalOnAll(poolobj, @clearvars, 0);
    end
    
    % Use the current 
    path = Folders{i_dir};
    
    % Optional: Prompt user for minimum central fitting timescale and the
    % indenter geometry
    % minTimescale = input('Please enter the minimum timescale to use for fitting (e.g. 1e-4): ');
    % tipOpts = {"spherical","conical"};
    % [indx,~] = listdlg('PromptString','Indenter Geometry',...
    %     'SelectionMode','single',...
    %     'ListString',tipOpts);
    % tipGeom = tipOpts{indx};
    
    % Alternatively, set those values manually
    tipGeom = "conical";                    % The experiment tip geometry for the files that are being loaded
    
    % Settings for how to use the loaded data during analysis
    useSmoothData = 0;                      % The user can choose to use filtered data, or the original raw data for fitting
    
    % Settings for loading the data
    loadDataSettings  = struct();
    
    % Required Settings:
    loadDataSettings.includeRetract = 0;             % Don't include data from the retract curve
    loadDataSettings.filterType = 'FIR';             % Choose the filter used to smooth data
    loadDataSettings.findRep = 'legendre';           % Search direction for the repulsive region
    loadDataSettings.removeNegatives = true;         % Remove negative values in the data stream
    loadDataSettings.createAverage = false;          % Create averaged rows at the END of the datastruct
    
    % Conditional Settings (depending on filter):
    loadDataSettings.N = 2;                          % Order of Butterworth filter (if used)
    loadDataSettings.cutoff_Hz = 5000;               % Cutoff frequency of Butterworth (if used)
    
    % Load the AFM Data
    dataStruct = LoadAFMData(path,loadDataSettings);
    
    % Allow for the contingency where multiple maps have been loaded from a
    % single folder by checking the number of unique values in mapID.
    numMaps = numel(unique([dataStruct(:).mapID]));
    
    if numMaps > 1
        % Grab the map filenames to name our output files easily
        FilesCheck=dir([originalPath filesep '*.*']);

        % Remove Directories
        FilesCheck=FilesCheck(~ismember({FilesCheck.name},{'.','..'}));
        toRemove = find([FilesCheck.isdir] == 1);
        FilesCheck(toRemove) = [];

        % Remove Filetypes To Ignore
        toRemove = find(~endsWith({FilesCheck.name}, {'.ibw','.txt','.spm','.mat','.csv','.jpk-qi-data'}));
        FilesCheck(toRemove) = [];

        toRemove = find(contains({FilesCheck.name}, {'settingsStruct','Settings'}));
        FilesCheck(toRemove) = [];

        toRemove = find(contains({FilesCheck.name}, {'FitResults','PlotResults','log.txt','mapStruct'}));
        FilesCheck(toRemove) = [];
        
        fileID = cell(length(FilesCheck),1);
        for i_file = 1:length(FilesCheck)
            temp = FilesCheck(i_file).name;
            temp = strsplit(temp, {'.'},'CollapseDelimiters',true);
            temp_name = temp{1};
            temp_name = strrep(temp_name,' ','');
            temp_name = strrep(temp_name,savePrependIn,'');
            fileID{i_file} = ['_' temp_name];
        end
    else
        % Fake identifier because which map this corresponds to is obvious
        fileID = {''};
    end
    
    for mapidx = 1:numMaps
        
        % Start Timer
        tic
        
        % Choose the right data from our structure
        if mapidx == 1
            ishift = 1;
            numPixels = dataStruct(ishift).mapSize(1)*dataStruct(ishift).mapSize(2);
            ai = 1;
            bi = numPixels;
        else
            ishift = ishift+numPixels;
            numPixels = dataStruct(ishift).mapSize(1)*dataStruct(ishift).mapSize(2);
            ai = ishift;
            bi = mapidx*numPixels;
        end

        % Create our cell arrays containing the data we care to send to the
        % ViscoFit class
        forces = cell(1,numPixels);
        times = cell(size(forces));
        indentations = cell(size(forces));
        tipSize = cell(size(forces));
        nu = cell(size(forces));
        mapSize = cell(size(forces));
        pixelHeight = cell(size(forces));

        ipix = 0;
        for i = ai:bi
            ipix = ipix+1;
            if loadDataSettings.removeNegatives
                if useSmoothData
                    forces{ipix} = dataStruct(i).F_r_smooth(dataStruct(i).F_r_smooth>0);
                    times{ipix} = dataStruct(i).t_r_smooth(dataStruct(i).F_r_smooth>0);
                    indentations{ipix} = dataStruct(i).h_r_smooth(dataStruct(i).F_r_smooth>0);
                else
                    forces{ipix} = dataStruct(i).F_r(dataStruct(i).F_r>0);
                    times{ipix} = dataStruct(i).t_r(dataStruct(i).F_r>0);
                    indentations{ipix} = dataStruct(i).h_r(dataStruct(i).F_r>0);
                end
            else
                if useSmoothData
                    forces{ipix} = dataStruct(i).F_r_smooth;
                    times{ipix} = dataStruct(i).t_r_smooth;
                    indentations{ipix} = dataStruct(i).h_r_smooth;
                else
                    forces{ipix} = dataStruct(i).F_r;
                    times{ipix} = dataStruct(i).t_r;
                    indentations{ipix} = dataStruct(i).h_r;
                end
            end
            tipSize{ipix} = dataStruct(i).r_tip;
            nu{ipix} = dataStruct(i).nu_sample;
            mapSize{ipix} = dataStruct(i).mapSize;
            pixelHeight{ipix} = dataStruct(i).height;
        end

        % Test the Fitting Functions using the Nelder-Mead Solver
        % Make a structure for our settings
        classSettings = struct;

        timescaleEst = 10^mode(ceil(log10(gradient([times{:}]))));

        % Test the Maxwell
        classSettings.minTimescale = 10*timescaleEst;       % Timescale for the first viscoelastic element
        classSettings.nu = nu;                              % Sample Poisson Ratio for all curves
        classSettings.tipGeom = tipGeom;                    % Tip geometry for these experiments
        classSettings.fitLog = false;                       % Log-scale-resample the data before fitting (faster)
        classSettings.thinSample = true;                    % Use the thin-sample correction feature
        classSettings.correctTilt = correctTilt;            % Correct any tilt in the substrate by fitting a plane to the image edges. Do NOT use for monolayers.
        classSettings.hideSubstrate = hideSubstrate;        % Remove pixels designated as "substrate" from the visco fitting
        
        if classSettings.correctTilt
            zeroSubstrate = true;                           % In addition to correcting for tilt, also set the new, flat surface to have a minimum value starting at zero.
            pixelHeight = fixMapTilt(mapSize,pixelHeight,zeroSubstrate);
        end
        
        classSettings.pixelHeight = pixelHeight;            % The height for each pixel in the map. Used for thresholding/thin sample correction

        % Create the class object
        viscoZ = ViscoFitZ(forces,times,indentations,tipSize,classSettings);

        % Make a structure for our settings
        fitSettings = struct;

        % Test the Maxwell
        fitSettings.solver = 'nelder-mead';                 % Fit using Nelder-Mead Simplex
        fitSettings.model = 'maxwell';                      % Use Generalized Maxwell Model
        fitSettings.n_elements = 3;                         % Fit iteratively for up to 3 elements
        fitSettings.elasticSetting = 1;                     % Include Elastic Term
        fitSettings.fluidSetting = 0;                       % No Steady-State Fluidity
        fitSettings.n_iterations = 1;                       % Use 10 random initializations
        fitSettings.n_fitIterations = 1e3;                  % No. of iterations for solver
        fitSettings.errortype = 'mse';                      % Use Mean-Squared Error during fitting
        fitSettings.N_workers = N_workers;                  % Pass the number of workers to the fitting
        fitSettings.hideSubstrate = viscoZ.hideSubstrate;   % Remove pixels designated as "substrate" from the visco fitting
        
        maxwellFit_zTransform = fitMapZ_func(viscoZ,fitSettings);
        maxwellFit_zTransform.mapSize = mode(cat(1,mapSize{:}),1);
        save([path filesep savePrepend fileID{mapidx} '_FitResults_NelderMead_zTransform.mat'],'maxwellFit_zTransform','-v7.3')
        clearvars maxwellFit_zTransform

%         % Test the Voigt
%         fitSettings.model = 'voigt';
%         voigtFit_zTransform = fitMapZ_func(visco,fitSettings);
%         voigtFit_zTransform.mapSize = mode(cat(1,mapSize{:}),1);
%         save([path filesep savePrepend fileID{mapidx} '_FitResults_NelderMead_zTransform.mat'],'voigtFit_zTransform','-append')
%         clearvars voigtFit_zTransform

%         % Test the PLR
%         fitSettings.model = 'PLR';
%         PLRFit_zTransform = fitMapZ_func(visco,fitSettings);
%         PLRFit_zTransform.mapSize = mode(cat(1,mapSize{:}),1);
%         save([path filesep savePrepend fileID{mapidx} '_FitResults_NelderMead_zTransform.mat'],'PLRFit_zTransform','-append')
%         clearvars PLRFit_zTransform
        
        clearvars forces times indentations tipSize nu mapSize pixelHeight
        
        analysisTime = toc;
        fprintf('\nTotal Map Analysis Time: %4.2f Minutes\n',analysisTime/60);
        
    end % End mapID for loop

    % Memory Management
    clearvars dataStruct
    
end % End folder for loop
    
%% Delete Parallel Pool of MATLAB Workers
poolobj = gcp('nocreate');
delete(poolobj);

end