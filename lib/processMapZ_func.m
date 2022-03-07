function fitStructOut = processMapZ_func(fitClass,varargin)
% PROCESSMAPZ_FUNC Extract Viscoelastic Information using the Z-Transform
%   This function takes in a variety of settings in addition to
%   the data already provided to the class and analyzes it. This specific
%   implementation utilizes the z-transform method.
%
%   This function will result in the observed viscoelastic information
%   for each pixel of a force map. To do this, the pixel
%   treatment is parallelized such that each worker handles the
%   extraction procedure for a group of pixels.

% Initialize Output Structure
tempFile = fullfile(tempdir,'fitStructTemp.mat');
if exist(tempFile,'file') == 2
    delete(tempFile);
end
fitStruct = matfile(tempFile,'Writable',true);

% Default Settings
N_workers = [];             % Number of workers for parpool
hideSubstrate = false;      % Setting to threshold based on height, and 
                            % remove the substrate from viscoelastic fitting. 
                            % Note: Do NOT use ignoreGlass on monolayers!
smoothOpt = 'none';         % Which smoothing setting to use on the harmonics.
                            % Options: none, g-time, ma-time, g-hz, ma-hz
                            % *-time smooths before z-transform. *-hz will
                            % smooth after z-transforming F and h.
windowsize = 0.05;          % Window size for smoothing methods

if ~isempty(varargin)
    % Only one varargin is accepted, and it is a structure
    % containing all of the settings information we require
    fitOpts = varargin{1};

    % Get the settings provided
    if isfield(fitOpts,'N_workers')
        N_workers = fitOpts.N_workers;
    end
    if isfield(fitOpts,'hideSubstrate')
        hideSubstrate = fitOpts.hideSubstrate;
    end
    if isfield(fitOpts,'smoothOpt')
        smoothOpt = fitOpts.smoothOpt;
    end
    if isfield(fitOpts,'windowsize')
        windowsize = fitOpts.windowsize;
    end
end

% Store the fit settings for future reference
fitStruct.hideSubstrate = hideSubstrate;
fitStruct.ViscoClass = fitClass;

% For Matlab's Parfor, we have to explicitly define the loop
% bounds ahead of time:
n_pixels = numel(fitClass.forces_cell);

% We initialize the map variables that we will be updating to have cell
% arrays of the appropriate size.
relaxanceMap = cell(size(fitClass.forces_cell));
retardanceMap = cell(size(fitClass.forces_cell));
alphaMap = cell(size(fitClass.forces_cell));
frequencyMap = cell(size(fitClass.forces_cell));
forceMap = cell(size(fitClass.forces_cell));
indMap = cell(size(fitClass.forces_cell));

% Create placeholders for the data we will obtain from our
% optimization attempts. These need to be pre-allocated to
% allow proper parallelization of the pixels.
fitStruct.relaxanceMap = relaxanceMap;
fitStruct.retardanceMap = retardanceMap;
fitStruct.alphaMap = alphaMap;
fitStruct.frequencyMap = frequencyMap;
fitStruct.forceMap = forceMap;
fitStruct.indMap = indMap;

% Additional Pre-allocating step to try and avoid growing arrays
% internally and preventin broadcasting of big variables to the main loop.
dataIn = cell(1,n_pixels);
for i = 1:n_pixels
    % Pre-Allocate
    relaxanceMap{i} = NaN(size(fitClass.forces_cell{i}));
    retardanceMap{i} = NaN(size(fitClass.forces_cell{i}));
    alphaMap{i} = NaN;
    frequencyMap{i} = NaN(size(fitClass.forces_cell{i}));
    
    % Make array input for parfor slicing (lowers memory usage because only
    % the i-th index is sent to a worker for the i-th iteration, and then
    % it is cleared before the next pixel is handled).
    dataIn{i} = {fitClass.times_cell{i},fitClass.dts_cell{i},fitClass.forces_cell{i},...
        fitClass.indentations_cell{i},fitClass.tipSize_cell{i},fitClass.nu_cell{i},...
        fitClass.tipGeom,{},fitClass.thinSample,fitClass.pixelHeight_cell{i}};
end

% Determine which pixels to ignore, if hideSubstrate == true
[minHeight,~] = min([fitClass.pixelHeight_cell{:}]);
substrateCutoff = minHeight + 100e-9;
pixelHeightArray = ([fitClass.pixelHeight_cell{:}]);
pixelOrder = 1:numel(fitClass.pixelHeight_cell);
pixelsToRemove = false(size(pixelHeightArray));
pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

pixelSkip = 1:numel([fitClass.pixelHeight_cell{:}]);
pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list

if hideSubstrate
    fprintf('\n%d Pixels of %d have been marked as "substrate" and will be skipped during analysis.\n', numel(pixelSkip), numel([fitClass.pixelHeight_cell{:}]))
end
    
% % Look at the data in 2D
% figure
% scatter(pixelOrder(~pixelsToRemove),pixelHeightArray(~pixelsToRemove),'bo');
% hold on
% scatter(pixelOrder(pixelsToRemove),pixelHeightArray(pixelsToRemove),'rx');
% hold off
% 
% [X,Y] = meshgrid(1:128,flip(1:128));
% pixelsToRemove3D = reshape(pixelsToRemove,128,128);
% Z = reshape(pixelHeightArray,128,128);
% 
% % See the REAL map data in 3D
% figure
% scatter3(X(~pixelsToRemove3D),Y(~pixelsToRemove3D),Z(~pixelsToRemove3D),10,'b')
% hold on
% scatter3(X(pixelsToRemove3D),Y(pixelsToRemove3D),Z(pixelsToRemove3D),10,'r')
% hold off

% Create parallel pool
poolFlag = false;           % Indicates that we have made a parpool here, and should clean it up after use.
if isempty(gcp('nocreate'))
    if isempty(N_workers)
        poolobj = parpool('IdleTimeout', 120);
    else
        poolobj = parpool(N_workers,'IdleTimeout', 120);
    end
    poolFlag = true;
end

% Don't perform fitting, only smoothing the Z-Transform results
fprintf('\nBeginning Map Analysis (No Fitting):\n');
parfor j = 1:n_pixels
    
    % Check if the pixel should be skipped
    if hideSubstrate
        if any(ismember(j,pixelSkip))
            % This is the "substrate" trigger! Skip this pixel
            % and enter in NaN for the output data.
            relaxanceMap(j) = {NaN};
            retardanceMap(j) = {NaN};
            alphaMap(j) = {NaN};
            frequencyMap(j) = {NaN};
            continue;
        end
    end
    
    % Check if this pixel is valid
    if any(isnan(dataIn{j}{3})) || isempty(dataIn{j}{3})
        % This is the "bad pixel" trigger! Skip this pixel
        % and enter in NaN for the output data.
        relaxanceMap(j) = {NaN};
        retardanceMap(j) = {NaN};
        alphaMap(j) = {NaN};
        frequencyMap(j) = {NaN};
        continue;
    end
    
    % Pass along the correct thinSample setting
    if any(ismember(j,pixelSkip))
        % This is the "substrate" trigger! If we are in this statement, the
        % pixel IS a substrate curve and we should NOT allow the
        % thinSample setting to be passed along.
        thinPixelLoop = false;
    else
        % This is the "cell surface" trigger! We will allow the user's
        % desired setting to be passed along.
        thinPixelLoop = dataIn{j}{9};
    end
    
    % Get Z-Transform Data from Curve
    [Q_hz,~,F_hz,~,h_hz,f_hz,alphaInit] = zTransformCurve(dataIn{j},smoothOpt,windowsize,thinPixelLoop);
    
    U_hz = 1./(Q_hz);
    relaxanceMap(j) = {Q_hz};
    retardanceMap(j) = {U_hz};
    alphaMap(j) = {alphaInit};
    frequencyMap(j) = {f_hz};
    forceMap(j) = {F_hz};
    indMap(j) = {h_hz};

end

fitStruct.relaxanceMap = relaxanceMap;
fitStruct.retardanceMap = retardanceMap;
fitStruct.alphaMap = alphaMap;
fitStruct.frequencyMap = frequencyMap;
fitStruct.forceMap = forceMap;
fitStruct.indMap = indMap;

fprintf('\nMap Analysis Complete!\n');

% Clean out the old parpool
if poolFlag
   % Get the current pool
    poolobj = gcp('nocreate');
    delete(poolobj); 
end

% Create output and clean up the temporary file
fitStructOut = load(fitStruct.Properties.Source,'-mat');
clearvars fitStruct
if exist(tempFile,'file') == 2
    delete(tempFile);
end

end

