classdef ViscoFitZ
    %ViscoFitZ Class Containing All Info Necessary to Extract Viscoelastic
    %Information using the Z-Transform Method
    %   This class, once initialized with the appropriate information, can
    %   perform a variety of viscoelastic parameterization operations. It
    %   must be initialized with the observed force, time, indentation, tip
    %   size (either spherical tip radius or indenter cone angle), the
    %   minimum timescale to begin fitting with: 
    %   "forces" - {1xN} - Cell array containing repulsive force array from
    %   experiments 1 to N (if more than one experiment is considered). The
    %   units are NEWTONS. The arrays must use DOUBLE PRECISION.
    
    %   "times" - {1xN} - Cell array containing repulsive time array from
    %   experiments 1 to N (if more than one experiment is considered). The
    %   units are SECONDS. The arrays must use DOUBLE PRECISION.
    
    %   "indentations" - {1xN} - Cell array containing repulsive
    %   indentation observed for experiments 1 to N (if more than one 
    %   experiment is considered). The units are METERS. The arrays must 
    %   use DOUBLE PRECISION.
    
    %   "tipSize" - {1xN} - Cell array containing the characteristic tip
    %   size, either the tip radius (for spherical indentation, the
    %   default) or cone angle (for conical indentation). Note that for
    %   conical indentation, an additional argument must be passed to
    %   overwrite the default spherical indentation setting. The units are
    %   METERS (spherical) or DEGREES (conical).  The values must use 
    %   DOUBLE PRECISION.
    
    %   Additionally, the class can take a variety of other optional
    %   settings:
    
    %   nu - {1xN} - Cell array of the sample's poisson's ratio, if known, 
    %   and not equal to 0.5 (incompressible) which is the default.
    
    %   pixelHeight - {1xN} - Cell array of the height for each curve, 
    %   which is used to calculate the bottom effect correction for soft
    %   samples on a hard substrate. This is NOT optional if the
    %   "thinSample" setting has been set to true!
    
    %   tipGeom - string containing the tip geometry, a string of either
    %   "spherical" or "conical", which determines how the tipSize argument
    %   is interpreted. The default is "spherical".
    
    %   thinSample - Logical (T/F) - the thin-sample correction setting, a
    %   logical (true/false), which enabled or disables the use of
    %   specialized LR_* functions that have additional height-thresholding
    %   and higher-order correction factors built into the data correction.
    
    %   correctTilt - Logical (T/F) - the tilt angle substrate correction 
    %   setting, logical (true/false), which enabled or disables the use of
    %   a linear regression plane-fit which will correct the substrate
    %   angle. Do NOT use this setting for monolayers where the substrate
    %   is not visible! It will remove parts of the map if this is the
    %   case.
    
    properties
        forces double
        forces_cell cell
        times double {mustBePositive}
        times_cell cell
        dts double {mustBePositive}
        dts_cell cell
        indentations double
        indentations_cell cell
        tipSize double {mustBePositive}
        tipSize_cell cell
        tipGeom string
        minTimescale double {mustBePositive}
        nu double {mustBePositive}
        nu_cell cell
        thinSample logical
        correctTilt logical
        hideSubstrate logical
        pixelHeight double
        pixelHeight_cell cell
    end
    
    methods (Static = false)
        % This region is for functions which reference the class (obj) in
        % their definition. This is the most convenient because it doesn't
        % require passing all of the data to each function call. However,
        % when the overhead is large with parallel processing, we do NOT
        % want the entire object as overhead to all workers because this
        % can seriously drain memory (see 128x128 maps!).
        
        function obj = ViscoFitZ(forces,times,indentations,tipSize,varargin)
            %ViscoFit Construct an instance of the ViscoFit class
            %   This class is utilized to perform the parameter
            %   optimization steps. It helps collect the relevant data for
            %   multiple load-level fitting (i.e. multiple approach
            %   velocities at the same point). The properties (forces,
            %   times, indentations, and tip size) must be supplied as
            %   cell arrays, where each entry corresponds to a different
            %   experiment. This is necessary for situations where the
            %   experiments are different lengths, to avoid needing to
            %   zero-pad the inputs and pass matrices to this class.
            
            %   Note, that unless specified in the varargin, the tipSize is
            %   assumed to be a tip radius (spherical indentation). If the
            %   second varargin is set to 'conical', then the contact
            %   mechanics are changed to account for the tipSize being the
            %   tip angle.

            if nargin >= 4
                
                % Handle optional arguments. Create default values, and modify
                % if the varargin contains the designated "real" values.
                tempnu = cell(size(forces));
                tempnu(:) = {0.5}; % Incompressible, 0.5
                temp = cellfun(@(nu,t) nu.*ones(size(t)),tempnu,times,'UniformOutput',false);
                obj.nu = horzcat(temp{:});
                obj.nu_cell = temp;
                if ~isempty(varargin)
                    if isa(varargin{1},'struct')
                        % Grab the settings structure
                        inputSettings = varargin{1};
                        
                        % Store the settings in the class
                        temp = cellfun(@(nu,t) nu.*ones(size(t)),inputSettings.nu,times,'UniformOutput',false);
                        obj.nu = horzcat(temp{:});
                        obj.nu_cell = temp;

                        if logical(inputSettings.thinSample)
                            temp = cellfun(@(h,t) h.*ones(size(t)),inputSettings.pixelHeight,times,'UniformOutput',false);
                            obj.pixelHeight = horzcat(temp{:});
                            obj.pixelHeight_cell = inputSettings.pixelHeight;
                        end
                        
                        if isfield(inputSettings,'tipGeom')
                            obj.tipGeom = string(inputSettings.tipGeom);
                        else
                            obj.tipGeom = "spherical";
                        end
                        
                        if isfield(inputSettings,'minTimescale')
                            obj.minTimescale = inputSettings.minTimescale;
                        else
                            obj.minTimescale = 1e-4;
                        end
                        
                        if isfield(inputSettings,'thinSample')
                            obj.thinSample = inputSettings.thinSample;
                        else
                            obj.thinSample = false;
                        end
                        
                        if isfield(inputSettings,'correctTilt')
                            obj.correctTilt = inputSettings.correctTilt;
                        else
                            obj.correctTilt = false;
                        end
                        
                        if isfield(inputSettings,'hideSubstrate')
                            obj.hideSubstrate = inputSettings.hideSubstrate;
                        else
                            obj.hideSubstrate = false;
                        end
                        
                    else
                        error('You are not passing the settings correctly to the ViscoFitZ Class Initialization. Please ensure the fifth argument is a structure containing your settings.');
                    end
                    
                end
                
                % Load the data streams
                obj.forces = horzcat(forces{:});
                obj.forces_cell = forces;

                obj.times = horzcat(times{:});
                obj.times_cell = times;

                obj.indentations = horzcat(indentations{:});
                obj.indentations_cell = indentations;

                try
                    temp = cellfun(@(t) round(mode(gradient(t)),1,'significant').*ones(size(t)),times,'UniformOutput',false);
                    obj.dts = horzcat(temp{:});
                    obj.dts_cell = temp;
                catch
                    disp('pause');
                end
                temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,times,'UniformOutput',false);
                obj.tipSize = horzcat(temp{:});
                obj.tipSize_cell = temp;
                
            else
                error('Not enough data provided to the class definition. Please supply at least the force, time, indentation, and radii.');
            end
            
        end
        
        function errorsum = SSE_Maxwell(obj,params,elasticSetting,fluidSetting,varargin)
            %SSE_Maxwell Calculate the SSE for the Maxwell model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Maxwell Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_Maxwell function. You can also pass an optional
            %   argument to this function, which is a string defining the
            %   error calculation method. The standard is 'sse' which does
            %   the sum of squared errors. Alternately, you can use the
            %   Mean Squared Error by specifying 'mse'. This will normalize
            %   each dataset individually by its variance and degrees of
            %   freedom.
            
            errortype = 'sse';
            if ~isempty(varargin)
                if length(varargin) <= 2
                    % Grab the error setting
                    errortype = varargin{1};
                else
                    error('You are passing too many arguments to SSE_maxwell. Please verify your inputs.');
                end
            end
            
            % Calculate test forces
            test_forces = LR_Maxwell(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom,elasticSetting,fluidSetting,obj.thinSample,obj.pixelHeight);
            
            % calculate global residual
            switch errortype
                case 'sse'
                    error_global = sum((obj.forces-test_forces).^2);
                    errorsum = sum(error_global);
                case 'mse'
                    normtemp = cellfun(@(ydata) (1./(movvar(ydata,3).^2))...
                        ./(length(ydata)-length(params)),obj.forces_cell,'UniformOutput',false);
                    normtemp = cell2mat(normtemp);
                    errorsum = sum(((test_forces-obj.forces).^2).*normtemp);
                otherwise
                    error('That error type is not implemented in SSE_maxwell. Please select an error type that exists in the function definition.')
            end
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e12;
            lb(modulusInds) = 1e-2;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;
            
            if length(params) > 1
                if fluidSetting
                    ub(2) = max(tauCenters)*1e2;
                    lb(2) = min(obj.dts);
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                errorsum = Inf;
            end
        end % End Maxwell SSE Function
        
        function [storageMod,lossMod,lossAngle] = harmonics_Maxwell(~,omega,fitStruct)
            %harmonics_Maxwell Calculate the Viscoelastic Harmonic
            %Quantities for a given Generalized Maxwell Model Parameter Set
                %   Calculate the Viscoelastic Storage Modulus, Loss
                %   Modulus, and Loss Angle (a.k.a., Loss Tangent) for the
                %   given set of Generalized Maxwell Model parameters
                %   acquired using the fitData() class function. The
                %   structure that is provided by the fitData() function
                %   must be included, in addition to the desired frequency
                %   array for evaluation.
                
            % Make sure we're using the right function
            if ~strcmpi(fitStruct.model,'maxwell')
                error('You attempted to pass parameters from a different viscoelastic model to the harmonics_Maxwell function. Please ensure you are passing the correct results.');
            end
                
            % Load the settings we need
            params = fitStruct.bestParams;
            elasticSetting = fitStruct.elasticSetting;
            fluidSetting = fitStruct.fluidSetting;
            
            % Make sure our input arrays are the right shape
            if isrow(params) params = params'; end
            if ~isrow(omega) omega = omega'; end
            
            % Make sure to erase any erroneous value contained within the
            % elastic parameter slot if there is no elastic term included.
            if ~elasticSetting
                params(1) = 0;
            end
            
            % Preallocate arrays/matrices
            storageMod = zeros(size(omega));
            lossMod = zeros(size(omega));
            stiffMatrix = repmat(params(3:2:end),1,length(omega));
            tauMatrix = repmat(params(4:2:end),1,length(omega));
            omegaMatrix = repmat(omega,length(params(3:2:end)),1);

            % Add the appropriate adjustment to the storage and loss moduli
            % according to the inclusion of steady-state fludity. When
            % fluidity is present, the elastic term must be included in the
            % loss modulus, and the storage modulus effects will no longer
            % be constant across all frequencies.
            if fluidSetting
                storageMod = storageMod + (params(1).*(omega.^2).*(params(2).^2))./(1+(omega.^2).*(params(2).^2));
                lossMod = lossMod + (params(1).*omega.*params(2))./(1+(omega.^2).*(params(2).^2));
            else
                storageMod = storageMod + params(1);
            end
            
            % Add the effect of the arms.
            storageMod = storageMod + sum((stiffMatrix.*(omegaMatrix.^2).*(tauMatrix.^2))./(1+(omegaMatrix.^2).*(tauMatrix.^2)),1);
            lossMod = lossMod + sum((stiffMatrix.*omegaMatrix.*tauMatrix)./(1+(omegaMatrix.^2).*(tauMatrix.^2)),1);
            lossAngle = atand(lossMod./storageMod);
                
        end % End of the Maxwell Harmonics Calculation Function
                
        function errorsum = SSE_Voigt(obj,params,elasticSetting,fluidSetting,varargin)
            %SSE_Voigt Calculate the SSE for the Voigt model
            %   Calculate the Sum of Squared Errors for the Generalized
            %   Voigt Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_Voigt function.
            
            errortype = 'sse';
            if ~isempty(varargin)
                if length(varargin) <= 2
                    % Grab the error setting
                    errortype = varargin{1};
                else
                    error('You are passing too many arguments to SSE_maxwell. Please verify your inputs.');
                end
            end
            
            % Calculate test indentation
            test_indentations = LR_Voigt(params,obj.times,obj.dts,obj.forces,obj.tipSize,obj.nu,obj.tipGeom,elasticSetting,fluidSetting,obj.thinSample,obj.pixelHeight);
            
            switch obj.tipGeom
                case "spherical"
                    beta = 1.5;
                case "conical"
                    beta = 2;
            end
            
            % calculate global residual
            switch errortype
                case 'sse'
                    sse_global = sum(((obj.indentations.^beta)-test_indentations).^2);
                    errorsum = sum(sse_global);
                case 'mse'
                    normtemp = cellfun(@(ydata) (1./(movvar(ydata,3).^2))...
                        ./(length(ydata)-length(params)),obj.indentations_cell,'UniformOutput',false);
                    normtemp = cell2mat(normtemp);
                    errorsum = sum(((test_indentations-(obj.indentations.^beta)).^2).*normtemp);
                otherwise
                    error('That error type is not implemented in SSE_voigt. Please select an error type that exists in the function definition.')
            end
            
            [tauInds,modulusInds] = getParamIndices(params);
            ub = zeros(size(params))+eps;
            lb = zeros(size(params));
            
            ub(modulusInds) = 1e2;
            lb(modulusInds) = 1e-12;
            
            tauCenters = obj.minTimescale.*(10.^( (1:length(params(3:2:end)))-1 ));
            ub(tauInds) = tauCenters*10;
            lb(tauInds) = tauCenters/10;

            if length(params) > 1
                if fluidSetting
                    ub(2) = 1;
                    lb(2) = 0;
                end
            end
            
            if any(ub-params < 0) || any(params-lb < 0)
                errorsum = Inf;
            end
        end % End Voigt SSE Map Function
        
        function [storageMod,lossMod,lossAngle] = harmonics_Voigt(~,omega,fitStruct)
            %harmonics_Voigt Calculate the Viscoelastic Harmonic
            %Quantities for a given Generalized Voigt Model Parameter Set
                %   Calculate the Viscoelastic Storage Modulus, Loss
                %   Modulus, and Loss Angle (a.k.a., Loss Tangent) for the
                %   given set of Generalized Voigt Model parameters
                %   acquired using the fitData() class function. The
                %   structure that is provided by the fitData() function
                %   must be included, in addition to the desired frequency
                %   array for evaluation.
            
            % Make sure we're using the right function
            if ~strcmpi(fitStruct.model,'voigt')
                error('You attempted to pass parameters from a different viscoelastic model to the harmonics_Voigt function. Please ensure you are passing the correct results.');
            end
                
            % Load the settings we need
            params = fitStruct.bestParams;
            elasticSetting = fitStruct.elasticSetting;
            fluidSetting = fitStruct.fluidSetting;
                
            % Make sure our input arrays are the right shape
            if isrow(params) params = params'; end
            if ~isrow(omega) omega = omega'; end
            
            % Make sure to erase any erroneous value contained within the
            % elastic parameter slot if there is no elastic term included.
            if ~elasticSetting
                params(1) = 0;
            end
            
            if ~fluidSetting
                params(2) = 0;
            end
            
            % Preallocate arrays/matrices
            storageCompliance = zeros(size(omega));
            lossCompliance = zeros(size(omega));
            compMatrix = repmat(params(3:2:end),1,length(omega));
            tauMatrix = repmat(params(4:2:end),1,length(omega));
            omegaMatrix = repmat(omega,length(params(3:2:end)),1);
            
            % Add the effects of the elastic arm and steady-state fluidity
            storageCompliance = storageCompliance + params(1);
            lossCompliance = lossCompliance + params(2)./omega;
            
            % Add the effect of the arms.
            storageCompliance = storageCompliance + sum(compMatrix./(1+(omegaMatrix.^2).*(tauMatrix.^2)),1);
            lossCompliance = lossCompliance + sum((compMatrix.*omegaMatrix.*tauMatrix)./(1+(omegaMatrix.^2).*(tauMatrix.^2)),1);

            % Calculate the absolute modulus
            absCompliance = sqrt(storageCompliance.^2 + lossCompliance.^2);
            
            % Add the effect of the arms.
            storageMod = storageCompliance./(absCompliance.^2);
            lossMod = lossCompliance./(absCompliance.^2);
            lossAngle = atand(lossMod./storageMod);

        end % End of the Voigt Harmonics Calculation Function
                
        function errorsum = SSE_PLR(obj,params,varargin)
            %SSE_PLR Calculate the SSE for the PLR model
            %   Calculate the Sum of Squared Errors for the Power Law
            %   Rheology Model according to the Lee and Radok indentation
            %   configuration, given a set of input parameters (params).
            %   This function performs fitting simultaneously for all force
            %   curves by linearizing the entire dataset (i.e. taking all
            %   curves for a particular approach velocity and stacking them
            %   end-to-end in a row vector). Separation and convolution of
            %   the correct regions of this row vector is handled inside
            %   the LR_PLR function.
            
            errortype = 'sse';
            erroropts = {'sse','mse'};
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if ischar(varargin{i}) && strcmpi(varargin{i},erroropts)
                        % Grab the error setting
                        errortype = varargin{i};
                    end
                end
            end
            
            % Calculate test forces
            test_forces = LR_PLR(params,obj.times,obj.dts,obj.indentations,obj.tipSize,obj.nu,obj.tipGeom,obj.thinSample,obj.pixelHeight);
            
            % calculate global residual
            switch errortype
                case 'sse'
                    sse_global = sum((obj.forces-test_forces).^2);
                    errorsum = sum(sse_global);
                case 'mse'
                    normtemp = cellfun(@(ydata) (1./(movvar(ydata,3).^2))...
                        ./(length(ydata)-length(params)),obj.forces_cell,'UniformOutput',false);
                    normtemp = cell2mat(normtemp);
                    errorsum = sum(((test_forces-obj.forces).^2).*normtemp);
                otherwise
                    error('That error type is not implemented in SSE_PLR. Please select an error type that exists in the function definition.')
            end
            
            % Power Law Rheology Roster:
            % [E_0 alpha]
            ub = [1e12;1];
            lb = [1e-2;0];
            
            if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
                errorsum = Inf;
            end
        end % End PLR SSE Function
        
        function [storageMod,lossMod,lossAngle] = harmonics_PLR(~,omega,fitStruct)
            %harmonics_PLR Calculate the Viscoelastic Harmonic
            %Quantities for a given Power Law Model Parameter Set
                %   Calculate the Viscoelastic Storage Modulus, Loss
                %   Modulus, and Loss Angle (a.k.a., Loss Tangent) for the
                %   given set of Power Law Rheology Model parameters
                %   acquired using the fitData() class function. The
                %   structure that is provided by the fitData() function
                %   must be included, in addition to the desired frequency
                %   array for evaluation. This function utilizes the DFT
                %   because the analytical Fourier transform of the PLR
                %   model is difficult to calculate. To provide a starting
                %   point for the manuscript, we have elected to estimate
                %   it discretely.
                
            % Make sure we're using the right function
            if ~strcmpi(fitStruct.model,'plr')
                error('You attempted to pass parameters from a different viscoelastic model to the harmonics_PLR function. Please ensure you are passing the correct results.');
            end
                
            % Load the settings we need
            params = fitStruct.bestParams;
            elasticSetting = fitStruct.elasticSetting;
            fluidSetting = fitStruct.fluidSetting;
            
            % Make sure our input arrays are the right shape
            if isrow(params) params = params'; end
            if ~isrow(omega) omega = omega'; end
            
            % Preallocate arrays/matrices
            storageMod = zeros(size(omega));
            lossMod = zeros(size(omega));
            
            % Calculate Power Law Rheology Modulus
            time_downsampled = round(fliplr(1./(omega./(2*pi))),4,'significant');
            time = time_downsampled(1):fitStruct.dt:time_downsampled(end);
            w_full = fliplr(2.*pi.*(1./time));
            t_prime = fitStruct.dt;
            E = params(1).*( (1 + time./t_prime) .^ (-params(2)) );
            G = E./(2*(1+fitStruct.nu_sample));
            
            % Calculate the discrete fourier transform (DFT) of the modulus
            % found using the PLR model.
            complexModulus = fft(G);
            storageMod = (real(complexModulus));
            storageMod(storageMod<1) = 1e-3;
            lossMod = (imag(complexModulus));
            lossMod(lossMod<1) = 1e-3;
            
            % Downsample again to the input frequency vector
            storageMod = interp1(w_full,storageMod,omega,'spline');
            lossMod = interp1(w_full,lossMod,omega,'spline');
            lossAngle = atand(lossMod./storageMod);
                
        end % End of the PLR Harmonics Calculation Function
        
        function [ub,lb] = getParamsCI(~,param_pop,varargin)
            %getParamsCI Get the Confidence Interval for a Parameter
            %Population
            %   Calculates the confidence bounds based upon the
            %   parameter population provided by the user in param_pop.
            %   If the user does not specify an accuracy through the
            %   varargin, the default is used (95% CI) using the
            %   student-t distribution (two-sided).
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    switch i
                        case 1
                            confInt = varargin{i};
                            if (confInt < 0) || (confInt > 1)
                                error('You provided an invalid accuracy level specification to getParamsCI(). Ensure your number remains between 0 and 1.')
                            end
                        otherwise
                            error('You have provided too many inputs to getParamsCI() through varargin.');
                    end
                end
            else
                confInt = 0.95;
            end
            
            % Quickly get our basic parameter population statistics
            meanPop = mean(param_pop,2);
            stdPop = std(param_pop,0,2);
            
            % Calculate the degrees of freedom in the system: N-p where N
            % is the number of observations (i.e., parameter sets in the
            % population) and p (the number of parameters in the model).
            dof = size(param_pop,2)-size(param_pop,1);
            
            % Get the critical values for our bounds. The way this is
            % written, the user can provide either the upper or lower bound
            % and still get the right scores. For example, if the user
            % gives 0.95 (for 95% CI's), the first element of tScores will
            % be bigger than the second (since the first has an input of
            % 0.95 and the second has 0.05). The values are then sorted in
            % increasing order, so the first element is smaller (used for
            % the lower bound). Alternately, if the user gives 0.05, then
            % the first element will be smaller than the second. The
            % absolute value is in place for just this case, as it will
            % therefore provide 0.05 and 0.95 to the two elements in
            % tScores.
            tScores = sort([tinv(abs(confInt),dof) tinv(abs(confInt-1),dof)]);
            
            % Calculate the upper and lower bounds, to be returned.
            lb = meanPop + tScores(1).*stdPop./sqrt(size(param_pop,2));
            ub = meanPop + tScores(2).*stdPop./sqrt(size(param_pop,2));
        end
        
        function fitStructOut = fitData(obj, varargin)
            % FITDATA Fit a Viscoelastic Model to Force Map Data with 
            %Z-Transform
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
            fitStructOut = matfile(tempFile,'Writable',true);

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
            fitStructOut.solver = solver;
            fitStructOut.model = model;
            fitStructOut.n_elements = n_elements;
            fitStructOut.elasticSetting = elasticSetting;
            fitStructOut.fluidSetting = fluidSetting;
            fitStructOut.n_iterations = n_iterations;
            fitStructOut.hideSubstrate = hideSubstrate;
            fitStructOut.ViscoClass = obj;

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
                    fitStructOut.fluidSetting = fluidSetting;
                case 'custom'
                    objFuncMap = @customFunc_Map_zTransform;
                otherwise
                    error('Your chosen solver-model combination is not implemented yet.');
            end

            % Create placeholders for the data we will obtain from our
            % optimization attempts. These need to be pre-allocated to
            % allow proper parallelization of the pixels.
            fitStructOut.bestParams = cell(1,n_elements);
            fitStructOut.paramPopulation = cell(1,n_elements);
            fitStructOut.paramPopulationResiduals = cell(1,n_elements);
            fitStructOut.upperParamCI = cell(1,n_elements);
            fitStructOut.lowerParamCI = cell(1,n_elements);
            fitStructOut.elasticFitTime = cell(1,n_elements);
            fitStructOut.fitTime = cell(1,n_elements);
            fitStructOut.relaxanceMap = cell(1,n_elements);
            fitStructOut.retardanceMap = cell(1,n_elements);
            fitStructOut.alphaMap = cell(1,n_elements);
            fitStructOut.frequencyMap = cell(1,n_elements);

            % For Matlab's Parfor, we have to explicitly define the loop
            % bounds ahead of time:
            n_pixels = numel(obj.forces_cell);

            % We initialize the map variables that we will be updating. For
            % the first iteration, these are all initially empty.
            % Afterward, they will contain the previous model
            % configuration's results. These are updated in turn and stored
            % before being overwritten.
            bestParamsMap = cell(size(obj.forces_cell));
            paramPopulationMap = cell(size(obj.forces_cell));
            paramPopulationResidualsMap = cell(size(obj.forces_cell));
            upperParamCIMap = cell(size(obj.forces_cell));
            lowerParamCIMap = cell(size(obj.forces_cell));
            elasticFitTimeMap = cell(size(obj.forces_cell));
            fitTimeMap = cell(size(obj.forces_cell));
            relaxanceMap = cell(size(obj.forces_cell));
            retardanceMap = cell(size(obj.forces_cell));
            alphaMap = cell(size(obj.forces_cell));
            frequencyMap = cell(size(obj.forces_cell));

            % Additional Pre-allocating step to try and avoid growing arrays and
            % creating our input datasets for each iteration of parfor.
            minTimescaleTemp = obj.minTimescale;
            lbFluidity = 10^( floor(min(log10(obj.dts)))+1 );
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
                relaxanceMap{i} = NaN(size(obj.forces_cell{i}));
                retardanceMap{i} = NaN(size(obj.forces_cell{i}));
                alphaMap{i} = NaN;
                frequencyMap{i} = NaN(size(obj.forces_cell{i}));

                % Input data
                dataIn{i} = {obj.times_cell{i},obj.dts_cell{i},obj.forces_cell{i},...
                    obj.indentations_cell{i},obj.tipSize_cell{i},obj.nu_cell{i},...
                    obj.tipGeom,obj.minTimescale,obj.thinSample,obj.pixelHeight_cell{i}};
            end

            % Determine which pixels to ignore, if hideSubstrate == true
            if hideSubstrate
                [minHeight,~] = min([obj.pixelHeight_cell{:}]);
                substrateCutoff = minHeight + 100e-9;
                pixelHeightArray = ([obj.pixelHeight_cell{:}]);
                pixelOrder = 1:numel(obj.pixelHeight_cell);
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel([obj.pixelHeight_cell{:}]);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list

                fprintf('\n%d Pixels of %d have been marked as "substrate" and will be skipped during fitting.\n', numel(pixelSkip), numel([obj.pixelHeight_cell{:}]))

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
                % function obj.getParamsCI without broadcasting the entire
                % class (with all the data in it!) to the parallel workers.
                % This minimizes overhead because this way the class is not
                % copied for each worker.
                for j = 1:n_pixels
                    [upperParamCIMap{j},lowerParamCIMap{j}] = obj.getParamsCI(bestParamsMap{j},0.95);
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
                    fitStructOut.relaxanceMap = relaxanceMap;
                    fitStructOut.retardanceMap = retardanceMap;
                    fitStructOut.alphaMap = alphaMap;
                    fitStructOut.frequencyMap = frequencyMap;
                end

                temp = fitStructOut.bestParams;
                temp{i} = bestParamsMap;
                fitStructOut.bestParams = temp;
                clearvars temp

                temp = fitStructOut.paramPopulation;
                temp{i} = paramPopulationMap;
                fitStructOut.paramPopulation = temp;
                clearvars temp

                temp = fitStructOut.paramPopulationResiduals;
                temp{i} = paramPopulationResidualsMap;
                fitStructOut.paramPopulationResiduals = temp;
                clearvars temp

                temp = fitStructOut.upperParamCI;
                temp{i} = upperParamCIMap;
                fitStructOut.upperParamCI = temp;
                clearvars temp

                temp = fitStructOut.lowerParamCI;
                temp{i} = lowerParamCIMap;
                fitStructOut.lowerParamCI = temp;
                clearvars temp

                temp = fitStructOut.elasticFitTime;
                temp{i} = elasticFitTimeMap;
                fitStructOut.elasticFitTime = temp;
                clearvars temp

                temp = fitStructOut.fitTime;
                temp{i} = fitTimeMap;
                fitStructOut.fitTime = temp;
                clearvars temp

                tmpSize = dir(tempFile);
                fprintf('complete. Current file size: %1.4g GB\n\n',tmpSize.bytes/(1e9));
                clearvars tmpSize

            end % End Iterative Term Introduction Loop

            clearvars bestParamsMap paramPopulationMap paramPopulationResidualsMap ...
                upperParamCIMap lowerParamCIMap elasticFitTimeMap ...
                fitTimeMap relaxanceMap retardanceMap alphaMap frequencyMap

            % Create output and clean up the temporary file
            fitStructOut = load(fitStructOut.Properties.Source,'-mat');
            clearvars fitStruct
            if exist(tempFile,'file')==2
                delete(tempFile);
            end
                        
        end % End fitData()
        
        function fitStructOut = processData(obj, varargin)
            % PROCESSDATA Process a QI Map using the Z-Transform Method 
            %   This function takes in a variety of settings in addition to
            %   the data already provided to the class and analyzes it. 
            %   This specific implementation utilizes the z-transform method.
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
            fitStruct.ViscoClass = obj;

            % For Matlab's Parfor, we have to explicitly define the loop
            % bounds ahead of time:
            n_pixels = numel(obj.forces_cell);

            % We initialize the map variables that we will be updating to have cell
            % arrays of the appropriate size.
            relaxanceMap = cell(size(obj.forces_cell));
            retardanceMap = cell(size(obj.forces_cell));
            alphaMap = cell(size(obj.forces_cell));
            frequencyMap = cell(size(obj.forces_cell));

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
                relaxanceMap{i} = NaN(size(obj.forces_cell{i}));
                retardanceMap{i} = NaN(size(obj.forces_cell{i}));
                alphaMap{i} = NaN;
                frequencyMap{i} = NaN(size(obj.forces_cell{i}));

                % Make array input for parfor slicing (lowers memory usage because only
                % the i-th index is sent to a worker for the i-th iteration, and then
                % it is cleared before the next pixel is handled).
                dataIn{i} = {obj.times_cell{i},obj.dts_cell{i},obj.forces_cell{i},...
                    obj.indentations_cell{i},obj.tipSize_cell{i},obj.nu_cell{i},...
                    obj.tipGeom,{},obj.thinSample,obj.pixelHeight_cell{i}};
            end

            % Determine which pixels to ignore, if hideSubstrate == true
            if hideSubstrate
                [minHeight,~] = min([obj.pixelHeight_cell{:}]);
                substrateCutoff = minHeight + 100e-9;
                pixelHeightArray = ([obj.pixelHeight_cell{:}]);
                pixelOrder = 1:numel(obj.pixelHeight_cell);
                pixelsToRemove = false(size(pixelHeightArray));
                pixelsToRemove(pixelHeightArray <= substrateCutoff) = true;

                pixelSkip = 1:numel([obj.pixelHeight_cell{:}]);
                pixelSkip(~pixelsToRemove) = [];    % Remove the pixels we want to keep from the list

                fprintf('\n%d Pixels of %d have been marked as "substrate" and will be skipped during fitting.\n', numel(pixelSkip), numel([obj.pixelHeight_cell{:}]))

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
            
        end % End processData()
        
    end % End Methods
    
    methods (Static = true)
        
        % Empty
        
    end % End Static Methods
    
end

