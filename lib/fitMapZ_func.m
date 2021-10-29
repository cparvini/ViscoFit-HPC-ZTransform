function fitStructOut = fitMapZ_func(fitClass,varargin)
% FITMAPZ_FUNC Fit a Viscoelastic Model to Force Map Data with Z-Transform
%   This function takes in a variety of settings in addition to
%   the data already provided to the class and performs an
%   optimization procedure based on those settings. This specific
%   implementation will utilize the z-transform method 
%
%   This function will result in an optimized parameter set
%   for each pixel of a force map. To do this, the pixel
%   treatment is parallelized such that each worker handles the
%   fitting procedure for a single force curve. As opposed to
%   having a single set of parameters for each model
%   configuration, as occurs for fitData(), each model
%   configuration will output a grid of parameter sets.
%   
%   In particular, the number of initializations chosen
%   (n_iterations) may need to be increased, as with the number
%   of solver iterations (n_fitIterations) and maximum number
%   of elements in the viscoelastic series (n_elements). This
%   function will iteratively introduce viscoelastic elements,
%   feeding forward the results from the previous iteration,
%   until the optimization has been performed for a model
%   configuration with a number of viscoelastic elements 
%   equal to n_elements.

% Initialize Output Structure
tempFile = fullfile(tempdir,'fitStructTemp.mat');
fitStruct = matfile(tempFile,'Writable',true);

% Default Settings
solver = 'nelder-mead';     % Fit using Nelder-Mead Simplex
model = 'maxwell';          % Use Generalized Maxwell Model
n_elements = 3;             % Fit iteratively for up to 3 elements
elasticSetting = 1;         % Include Elastic Term
fluidSetting = 0;           % No Steady-State Fluidity
n_iterations = 10;          % Use 10 random initializations
n_fitIterations = 1e4;      % No. of iterations for solver
N_workers = [];             % Number of workers for parpool
hideSubstrate = false;      % Setting to threshold based on height, and 
                            % remove the substrate from viscoelastic fitting. 
                            % Note: Do NOT use ignoreGlass on monolayers!
smoothOpt = 'ma-time';      % Which smoothing setting to use on the harmonics.
                            % Options: none, g-time, ma-time, g-hz, ma-hz
                            % *-time smooths before z-transform. *-hz will
                            % smooth after z-transforming F and h.
windowsize = 10;            % Window size for smoothing methods
if ~isempty(varargin)
    % Only one varargin is accepted, and it is a structure
    % containing all of the settings information we require
    fitOpts = varargin{1};

    % Get the settings provided
    try
        solver = lower(fitOpts.solver);
    end
    try
        model = lower(fitOpts.model);
    end
    try
        n_elements = fitOpts.n_elements;
    end
    try
        elasticSetting = fitOpts.elasticSetting;
    end
    try
        fluidSetting = fitOpts.fluidSetting;
    end
    try
        n_iterations = fitOpts.n_iterations;
    end
    try
        n_fitIterations = fitOpts.n_fitIterations;
    end
    try
        N_workers = fitOpts.N_workers;
    end
    try
        hideSubstrate = fitOpts.hideSubstrate;
    end
    if strcmpi(solver,'custom')
        if isfield(fitOpts,'customFunc')
            customFunc = fitOpts.customFunc;
        else
            error('You selected "custom" for the solver type and did not supply the function name in the class settings.');
        end
    end
end

% Store the fit settings for future reference
fitStruct.solver = solver;
fitStruct.model = model;
fitStruct.n_elements = n_elements;
fitStruct.elasticSetting = elasticSetting;
fitStruct.fluidSetting = fluidSetting;
fitStruct.n_iterations = n_iterations;
fitStruct.hideSubstrate = hideSubstrate;
fitStruct.ViscoClass = fitClass;

% Get the correct objective function for optimization
% Note we are referencing external functions here, so we use
% can avoid broadcasting the entire class just to run our
% objective map function! Each function has been defined in the
% class as a static method, too, for future reference but it is
% NOT recommended to use them for the purpose of OPTIMIZATION!
switch model
    case 'maxwell'
        objFuncMap = @SSE_Maxwell_Map_zTransform;
    case 'voigt'
        objFuncMap = @SSE_Voigt_Map_zTransform;
    case 'plr'
        objFuncMap = @SSE_PLR_Map_zTransform;
        % The storage roster for plr is different, and requires
        % the second index to maintain consistency. Thus, we
        % have to manually force the second term to be fit by
        % including the fluidity. This second position actually
        % corresponds to the exponent, alpha.
        fluidSetting = 1;
        fitStruct.fluidSetting = fluidSetting;
    case 'custom'
        objFuncMap = @customFunc_Map_zTransform;
    otherwise
        error('Your chosen solver-model combination is not implemented yet.');
end

% Create placeholders for the data we will obtain from our
% optimization attempts. These need to be pre-allocated to
% allow proper parallelization of the pixels.
fitStruct.bestParams = cell(1,n_elements);
fitStruct.paramPopulation = cell(1,n_elements);
fitStruct.paramPopulationResiduals = cell(1,n_elements);
fitStruct.upperParamCI = cell(1,n_elements);
fitStruct.lowerParamCI = cell(1,n_elements);
fitStruct.elasticFitTime = cell(1,n_elements);
fitStruct.fitTime = cell(1,n_elements);
fitStruct.relaxanceMap = cell(1,n_elements);
fitStruct.retardanceMap = cell(1,n_elements);
fitStruct.alphaMap = cell(1,n_elements);
fitStruct.frequencyMap = cell(1,n_elements);
        
% For Matlab's Parfor, we have to explicitly define the loop
% bounds ahead of time:
n_pixels = numel(fitClass.forces_cell);

% We initialize the map variables that we will be updating. For
% the first iteration, these are all initially empty.
% Afterward, they will contain the previous model
% configuration's results. These are updated in turn and stored
% before being overwritten.
bestParamsMap = cell(size(fitClass.forces_cell));
paramPopulationMap = cell(size(fitClass.forces_cell));
paramPopulationResidualsMap = cell(size(fitClass.forces_cell));
upperParamCIMap = cell(size(fitClass.forces_cell));
lowerParamCIMap = cell(size(fitClass.forces_cell));
elasticFitTimeMap = cell(size(fitClass.forces_cell));
fitTimeMap = cell(size(fitClass.forces_cell));
relaxanceMap = cell(size(fitClass.forces_cell));
retardanceMap = cell(size(fitClass.forces_cell));
alphaMap = cell(size(fitClass.forces_cell));
frequencyMap = cell(size(fitClass.forces_cell));

% Additional Pre-allocating step to try and avoid growing arrays and
% creating our input datasets for each iteration of parfor.
minTimescaleTemp = fitClass.minTimescale;
lbFluidity = 10^( floor(min(log10(fitClass.dts)))+1 );
dataIn = cell(1,n_pixels);
for i = 1:n_pixels
    % Pre-allocation
    bestParamsMap{i} = NaN(4,1);
    paramPopulationMap{i} = NaN(4,n_iterations);
    paramPopulationResidualsMap{i} = NaN(1,n_iterations);
    upperParamCIMap{i} = NaN(4,1);
    lowerParamCIMap{i} = NaN(4,1);
    elasticFitTimeMap{i} = NaN;
    fitTimeMap{i} = NaN;
    relaxanceMap{i} = NaN(size(fitClass.forces_cell{i}));
    retardanceMap{i} = NaN(size(fitClass.forces_cell{i}));
    alphaMap{i} = NaN;
    frequencyMap{i} = NaN(size(fitClass.forces_cell{i}));
    
    % Input data
    dataIn{i} = {fitClass.times_cell{i},fitClass.dts_cell{i},fitClass.forces_cell{i},...
        fitClass.indentations_cell{i},fitClass.tipSize_cell{i},fitClass.nu_cell{i},...
        fitClass.tipGeom,fitClass.minTimescale,fitClass.thinSample,fitClass.pixelHeight_cell{i}};
end

% Determine which pixels to ignore, if hideSubstrate == true
if hideSubstrate
    [minHeight,~] = min([fitClass.pixelHeight_cell{:}]);
    substrateCutoff = minHeight + 100e-9;
    pixelHeightArray = ([fitClass.pixelHeight_cell{:}]);
    pixelOrder = 1:numel(fitClass.pixelHeight_cell);
    pixelsToRemove = false(size(pixelHeightArray));
    pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;
    
    pixelSkip = 1:numel([fitClass.pixelHeight_cell{:}]);
    pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list
    
    fprintf('\n%d Pixels of %d have been marked as "substrate" and will be skipped during fitting.\n', numel(pixelSkip), numel([fitClass.pixelHeight_cell{:}]))
    
%     % Look at the data in 2D
%     figure
%     scatter(pixelOrder(~pixelsToRemove),pixelHeightArray(~pixelsToRemove),'bo');
%     hold on
%     scatter(pixelOrder(pixelsToRemove),pixelHeightArray(pixelsToRemove),'rx');
%     hold off
%     
%     [X,Y] = meshgrid(1:128,flip(1:128));
%     pixelsToRemove3D = reshape(pixelsToRemove,128,128);
%     Z = reshape(pixelHeightArray,128,128);
%     
%     % See the REAL map data in 3D
%     figure
%     scatter3(X(~pixelsToRemove3D),Y(~pixelsToRemove3D),Z(~pixelsToRemove3D),10,'b')
%     hold on
%     scatter3(X(pixelsToRemove3D),Y(pixelsToRemove3D),Z(pixelsToRemove3D),10,'r')
%     hold off

else
    
    % No pixels to skip!
    pixelSkip = [];
    
end

% Begin the iterative term introduction loop
fprintf('\nBeginning Map Analysis (%s):\n', model)
for i = 1:n_elements

    if ~isempty(gcp('nocreate'))
       % Get the current pool
        poolobj = gcp('nocreate');

        % Clean up the workers (memory management).
        % Unfortunately this can ONLY be done by restarting the
        % entire pool. This doesn't happen often for us, so
        % this is fine, but as soon as Mathworks will let users
        % clear the parpool of the junk from the previous
        % parfor loop, that fix should be used here instead
        % because starting a parallel pool is time intensive.
        delete(poolobj);

        % The only reason we have this here, as well as after
        % the loop, is because if there is a pool existing in
        % the main script we want to kill that before we start
        % the first iteration. This if statement will not
        % activate for parallel pools made in this for loop,
        % since the pool is deleted right after the loop
        % (memory management).

        % If the pool is NOT restarted, the memory will not be
        % liberated between parfor loops, causing MASSIVE
        % buildup. This means even a small amount of fitting
        % data (<1GB output file) can exceed 250 GB of
        % RAM/Memory during the first term search. Not good!
    end

    % Make a fresh pool
    if isempty(N_workers)
        poolobj = parpool('IdleTimeout', 120);
    else
        poolobj = parpool(N_workers,'IdleTimeout', 120);
    end

    % Send the class to the workers
    addAttachedFiles(poolobj, {'LR_Maxwell.m','LR_Voigt.m','LR_PLR.m',...
        'SSE_Maxwell_Map.m','SSE_Voigt_Map.m','SSE_PLR_Map.m',...
        'fitPixelZ.m'})

    fprintf('%d terms...', i)

    % Initialize Variables
    paramLen = length(bestParamsMap{1});
    
    % Original
    parfor j = 1:n_pixels
        % Check if the pixel should be skipped
        if hideSubstrate
            if any(ismember(j,pixelSkip))
                % This is the "ignore pixel" trigger! Skip this pixel
                % and enter in NaN for the output data.
                elasticFitTimeMap(j) = {0};
                fitTimeMap(j) = {0};
                if i > 1
                    bestParamsMap(j) = {NaN(paramLen+2,1)};
                    paramPopulationMap(j) = {NaN(paramLen+2,1)};
                else
                    bestParamsMap(j) = {NaN(paramLen,1)};
                    paramPopulationMap(j) = {NaN(paramLen,1)};
                end
                paramPopulationResidualsMap(j) = {NaN};
                if i == 1
                    relaxanceMap(j) = {NaN};
                    retardanceMap(j) = {NaN};
                    alphaMap(j) = {NaN};
                    frequencyMap(j) = {NaN};
                end
                continue;
            end
        end
        
        % Check if this pixel is valid
        if any(isnan(dataIn{j}{3})) || isempty(dataIn{j}{3})
            elasticFitTimeMap(j) = {0};
            fitTimeMap(j) = {0};
            if i > 1
                bestParamsMap(j) = {NaN(paramLen+2,1)};
                paramPopulationMap(j) = {NaN(paramLen+2,1)};
            else
                bestParamsMap(j) = {NaN(paramLen,1)};
                paramPopulationMap(j) = {NaN(paramLen,1)};
            end
            paramPopulationResidualsMap(j) = {NaN};
            if i == 1
                relaxanceMap(j) = {NaN};
                retardanceMap(j) = {NaN};
                alphaMap(j) = {NaN};
                frequencyMap(j) = {NaN};
            end
            continue;
        end
        
        [relaxanceMap(j),...
            retardanceMap(j),...
            alphaMap(j),...
            frequencyMap(j),...
            bestParamsMap(j),...
            paramPopulationMap(j),...
            paramPopulationResidualsMap(j),...
            elasticFitTimeMap(j),...
            fitTimeMap(j)] = fitPixelZ(i,n_iterations,model,solver,n_fitIterations,minTimescaleTemp,dataIn(j),paramLen,objFuncMap,elasticSetting,fluidSetting,bestParamsMap{j});
        
    end

    % Run through and calculate CI for parameters. This has to
    % be done here with another loop because we can't call the
    % function fitClass.getParamsCI without broadcasting the entire
    % class (with all the data in it!) to the parallel workers.
    % This minimizes overhead because this way the class is not
    % copied for each worker.
    for j = 1:n_pixels
        [upperParamCIMap{j},lowerParamCIMap{j}] = fitClass.getParamsCI(bestParamsMap{j},0.95);
    end

    if any(strcmp(model,'plr'))
        % The PLR model does not utilize more than a single
        % term, so "iterative term introduction" has no
        % meaning in this case. As such, the remaining
        % iterations are skipped and the results are returned

        % To stop the iterative term introduction for one of
        % the models, add an additional strcmp(model,'{name}')
        % after the one in this if statement declaration,
        % separated by a comma. The any() will trigger if the
        % type is PLR or whatever new entry you have added.
        % This is primarily in the case where you have a
        % custom model that does not use iterative term
        % introduction.
        break;

    end

    if ~isempty(gcp('nocreate'))
       % Get the current pool
        poolobj = gcp('nocreate');

        % Clean up the workers (memory management).
        % Unfortunately this can ONLY be done by restarting the
        % entire pool. This doesn't happen often for us, so
        % this is fine, but as soon as Mathworks will let users
        % clear the parpool of the junk from the previous
        % parfor loop, that fix should be used here instead
        % because starting a parallel pool is time intensive.
        delete(poolobj); 
    end

    % Store updated values in output structure
    if i == 1
        fitStruct.relaxanceMap = relaxanceMap;
        fitStruct.retardanceMap = retardanceMap;
        fitStruct.alphaMap = alphaMap;
        fitStruct.frequencyMap = frequencyMap;
    end
    
    temp = fitStruct.bestParams;
    temp{i} = bestParamsMap;
    fitStruct.bestParams = temp;
    clearvars temp

    temp = fitStruct.paramPopulation;
    temp{i} = paramPopulationMap;
    fitStruct.paramPopulation = temp;
    clearvars temp

    temp = fitStruct.paramPopulationResiduals;
    temp{i} = paramPopulationResidualsMap;
    fitStruct.paramPopulationResiduals = temp;
    clearvars temp

    temp = fitStruct.upperParamCI;
    temp{i} = upperParamCIMap;
    fitStruct.upperParamCI = temp;
    clearvars temp

    temp = fitStruct.lowerParamCI;
    temp{i} = lowerParamCIMap;
    fitStruct.lowerParamCI = temp;
    clearvars temp

    temp = fitStruct.elasticFitTime;
    temp{i} = elasticFitTimeMap;
    fitStruct.elasticFitTime = temp;
    clearvars temp

    temp = fitStruct.fitTime;
    temp{i} = fitTimeMap;
    fitStruct.fitTime = temp;
    clearvars temp

    tmpSize = dir(tempFile);
    fprintf('complete. Current file size: %1.4g GB\n\n',tmpSize.bytes/(1e9));
    clearvars tmpSize

end % End Iterative Term Introduction Loop

clearvars bestParamsMap paramPopulationMap paramPopulationResidualsMap ...
    upperParamCIMap lowerParamCIMap elasticFitTimeMap ...
    fitTimeMap relaxanceMap retardanceMap alphaMap frequencyMap

% Create output and clean up the temporary file
fitStructOut = load(fitStruct.Properties.Source,'-mat');
clearvars fitStruct
if exist(tempFile,'file')==2
    delete(tempFile);
end

end

