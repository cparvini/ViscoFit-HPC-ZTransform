function fitStructOut = processMapZ_func(fitClass,varargin)
% PROCESSMAPZ_FUNC Fit a Viscoelastic Model to Force Map Data with 
%Z-Transform
%   This function takes in a variety of settings in addition to
%   the data already provided to the class and analyzes it. This specific
%   implementation utilizes the z-transform method.
%
%   This function will result in the extracted viscoelastic info
%   for each pixel of a force map. To do this, the pixel
%   treatment is parallelized such that each worker handles the
%   exrtaction procedure for a single force curve.

% Initialize Output Structure
tempFile = fullfile(tempdir,'fitStructTemp.mat');
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
windowsize = 10;            % Window size for smoothing methods

if ~isempty(varargin)
    % Only one varargin is accepted, and it is a structure
    % containing all of the settings information we require
    fitOpts = varargin{1};

    % Get the settings provided
    try
        N_workers = fitOpts.N_workers;
    end
    try
        hideSubstrate = fitOpts.hideSubstrate;
    end
    try
        smoothOpt = fitOpts.smoothOpt;
    end
    try
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

% Create placeholders for the data we will obtain from our
% optimization attempts. These need to be pre-allocated to
% allow proper parallelization of the pixels.
fitStruct.relaxanceMap = relaxanceMap;
fitStruct.retardanceMap = retardanceMap;
fitStruct.alphaMap = alphaMap;
fitStruct.frequencyMap = frequencyMap;

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

    % Solve for z-transform quantities.
    % First, we need to find the correct alpha value and store
    % it for later. Create the initial guess first, using the
    % values from the force array.
    tipSize = dataIn{j}{5};
    nu = dataIn{j}{6};
    tipGeom = dataIn{j}{7};
    thinSample = dataIn{j}{9};
    h_finite = dataIn{j}{10};

    alphaInit = NaN(2,1);
    alphaInit(1) = (1/numel(dataIn{j}{1}))*(real(log(dataIn{j}{3}(1))) ...
        / real(log(dataIn{j}{3}(end))));
    alphaInit(2) = (1/numel(dataIn{j}{1}))*(real(log(dataIn{j}{4}(1))) ...
        / real(log(dataIn{j}{4}(end))));
    alphaInit = max(alphaInit(~isinf(alphaInit)),[],'omitnan');

    c = NaN;
    betaParam = NaN;
    switch tipGeom
        case "spherical"
            c = (8*sqrt(tipSize))./(3*(1-nu));
            betaParam = 1.5;
        case "conical"
            c = (2.*tan(tipSize.*pi./180))./(pi.*(1-nu.^2));
            betaParam = 2;
    end

    % Optimize alpha
    modelfun = @(alpha,x) (1/numel(x)).*fftshift(fft(x.*(exp(-alpha.*(1:numel(x))))));
    F_hz = NaN;
    h_hz = NaN;
    if ~thinSample
        F_t = dataIn{j}{3}; % forces
        h_t = dataIn{j}{4}.^(betaParam); % indentations
        
        switch smoothOpt
            case 'none'
                % Do not smooth
                F_hz = modelfun(alphaInit,F_t);
                h_hz = modelfun(alphaInit,h_t);
                
            case 'g-time'
                % Gaussian smoothing in time
                F_hz = modelfun(alphaInit,smoothdata(F_t,'gaussian',windowsize));
                h_hz = modelfun(alphaInit,smoothdata(h_t,'gaussian',windowsize));
                                
            case 'ma-time'
                % Moving Average smoothing in time
                F_hz = modelfun(alphaInit,smooth(dataIn{j}{1},F_t,windowsize)');
                h_hz = modelfun(alphaInit,smooth(dataIn{j}{1},h_t,windowsize)');
                                
            case 'g-hz'
                % Gaussian smoothing in frequency
                F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowsize);
                h_hz = smoothdata(modelfun(alphaInit,h_t),'gaussian',windowsize);
                                
            case 'ma-hz'
                % Moving Average smoothing in frequency
                F_hz = smooth(dataIn{j}{1},modelfun(alphaInit,F_t),windowsize)';
                h_hz = smooth(dataIn{j}{1},modelfun(alphaInit,h_t),windowsize)';                
                
        end
        
        dt = mode(dataIn{j}{2});
        fs = dt.^-1;
        N = length(F_hz);
        df = fs/N;
        f_hz = -fs/2:df:fs/2-df + (df/2)*mod(N,2);
        
        Q_hz = (F_hz./h_hz) .* (1./c);
        
%         % Observe Smoothing Results!
%         try
%             figure(hfig)
%             clf;
%         catch
%             hfig = figure;
%             tiledlayout(2,1);
%         end
%         nexttile
%         scatter(dataIn{j}{1},dataIn{j}{3}./BEC,'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},F_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{j}{1},dataIn{j}{4}.^(betaParam),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},h_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Indentation [m]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         try
%             figure(hzfig)
%             clf;
%         catch
%             hzfig = figure;
%             tiledlayout(2,1);
%         end
%         nexttile
%         scatter(dataIn{j}{1},modelfun(alphaInit,dataIn{j}{3}./BEC),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},F_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{j}{1},modelfun(alphaInit,dataIn{j}{4}.^(betaParam)),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},h_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Indentation [m]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         try
%             figure(Qhzfig)
%             clf;
%         catch
%             Qhzfig = figure;
%         end
%         scatter(f_hz,(modelfun(alphaInit,dataIn{j}{3}./BEC)./modelfun(alphaInit,dataIn{j}{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%         hold on
%         plot(f_hz,Q_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         xlim([min(f_hz) max(f_hz)])
%         grid on
%         hold off
%         legend('Original','Smooth','location','best')
        
    else
        
        h_t = dataIn{j}{4}.^(betaParam); % indentations
        BEC = NaN;
        
        switch tipGeom
            case "spherical"
                % Defined per Dimitriadis et al. (Biophy. Journ., 2002)
                % Assume sample adhered to surface
                chi_h = sqrt(tipSize.*h_t)./h_finite;
                alpha0 = -((1.2876 - 1.4678.*nu + 1.3442.*(nu.^2))./(1 - nu));
                beta0 = -((0.6387 - 1.0277.*nu + 1.5164.*(nu.^2))./(1 - nu));
                BEC = 1 - (2.*alpha0./pi).*chi_h ...
                    + (4.*(alpha0.^2)./(pi^2)).*(chi_h.^2) ...
                    - (8/(pi^3)).*((alpha0.^3)+(4*(pi^2)/15).*beta0).*(chi_h^3) ...
                    + (16.*alpha0./(pi^4)).*...
                    ((alpha0.^3)+(3*(pi^2)/5).*beta0).*(chi_h^4);
                
            case "conical"
                % Defined per Gavara & Chadwick (Nat. Nanotech. 2012)
                % Assume sample adhered to surface
                zeta = 1.7795;
                BEC = 1 + zeta.*((2.*tan(tipSize.*pi./180).*h_t)./((pi^2).*h_finite))...
                    + 16.*(zeta.^2).*(tan(tipSize.*pi./180).^2).*((h_t.^2)./(h_finite.^2))...
                    + 0; % Can include higher orders if we can find them somewhere.
        end
        
        F_t = dataIn{j}{3}./BEC; % forces
                        
        % Perform Smoothing
%         disp('pause');
        
        switch smoothOpt
            case 'none'
                % Do not smooth
                F_hz = modelfun(alphaInit,F_t);
                h_hz = modelfun(alphaInit,h_t);
                
            case 'g-time'
                % Gaussian smoothing in time
                F_hz = modelfun(alphaInit,smoothdata(F_t,'gaussian',windowsize));
                h_hz = modelfun(alphaInit,smoothdata(h_t,'gaussian',windowsize));
                                
            case 'ma-time'
                % Moving Average smoothing in time
                F_hz = modelfun(alphaInit,smooth(dataIn{j}{1},F_t,windowsize)');
                h_hz = modelfun(alphaInit,smooth(dataIn{j}{1},h_t,windowsize)');
                                
            case 'g-hz'
                % Gaussian smoothing in frequency
                F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowsize);
                h_hz = smoothdata(modelfun(alphaInit,h_t),'gaussian',windowsize);
                                
            case 'ma-hz'
                % Moving Average smoothing in frequency
                F_hz = smooth(dataIn{j}{1},modelfun(alphaInit,F_t),windowsize)';
                h_hz = smooth(dataIn{j}{1},modelfun(alphaInit,h_t),windowsize)';                
                
        end
        
        dt = mode(dataIn{j}{2});
        fs = dt.^-1;
        N = length(F_hz);
        df = fs/N;
        f_hz = -fs/2:df:fs/2-df + (df/2)*mod(N,2);
        
        Q_hz = (F_hz./h_hz) .* (1./c);
        
%         % Observe Smoothing Results!
%         try
%             figure(hfig)
%             clf;
%         catch
%             hfig = figure;
%             tiledlayout(2,1);
%         end
%         nexttile
%         scatter(dataIn{j}{1},dataIn{j}{3}./BEC,'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},F_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{j}{1},dataIn{j}{4}.^(betaParam),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},h_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Indentation [m]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         try
%             figure(hzfig)
%             clf;
%         catch
%             hzfig = figure;
%             tiledlayout(2,1);
%         end
%         nexttile
%         scatter(dataIn{j}{1},modelfun(alphaInit,dataIn{j}{3}./BEC),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},F_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{j}{1},modelfun(alphaInit,dataIn{j}{4}.^(betaParam)),'bo','linewidth',3)
%         hold on
%         plot(dataIn{j}{1},h_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Indentation [m]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         try
%             figure(Qhzfig)
%             clf;
%         catch
%             Qhzfig = figure;
%         end
%         scatter(f_hz,(modelfun(alphaInit,dataIn{j}{3}./BEC)./modelfun(alphaInit,dataIn{j}{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%         hold on
%         plot(f_hz,Q_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         xlim([min(f_hz) max(f_hz)])
%         grid on
%         hold off
%         legend('Original','Smooth','location','best')

    end

    U_hz = 1./(Q_hz);
    relaxanceMap(j) = {Q_hz};
    retardanceMap(j) = {U_hz};
    alphaMap(j) = {alphaInit};
    frequencyMap(j) = {f_hz};  

end

fitStruct.relaxanceMap = relaxanceMap;
fitStruct.retardanceMap = retardanceMap;
fitStruct.alphaMap = alphaMap;
fitStruct.frequencyMap = frequencyMap;

% Create output and clean up the temporary file
fitStructOut = load(fitStruct.Properties.Source,'-mat');
clearvars fitStruct
if exist(tempFile,'file')==2
    delete(tempFile);
end

end

