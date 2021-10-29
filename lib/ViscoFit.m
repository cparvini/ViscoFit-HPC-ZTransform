classdef ViscoFit
    %ViscoFit Class Containing All Info Necessary to Extract Viscoelastic
    %Parameters
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
    
    %   "minTimescale" - double - Single value which will be the "center"
    %   of the characteristic time range for the first viscoelastic element
    %   in the generalized rheological models. This is not necessary when
    %   using the PLR model; when using PLR, set this value to ONE.  The 
    %   value must use DOUBLE PRECISION.
    
    %   Additionally, the class can take a variety of other settings:
    
    %   Optional Argument #1 - the sample's poisson's ratio, if known, and
    %   not equal to 0.5 (incompressible) which is the default.
    
    %   Optional Argument #2 - the tip geometry, a string of either
    %   "spherical" or "conical", which determines how the tipSize argument
    %   is interpreted.
    
    %   Optional Argument #3 - the thin-sample correction setting, a
    %   logical (true/false), which enabled or disables the use of
    %   specialized LR_* functions that have additional height-thresholding
    %   and higher-order correction factors built into the data correction.
    
    %   Optional Argument #4 - the tilt angle substrate correction setting,
    %   logical (true/false), which enabled or disables the use of
    %   a linear regression plane-fit which will correct the substrate
    %   angle. Do NOT use this setting for monolayers where the substrate
    %   is not visible!
    
    %   Optional Argument #5 - the thin-sample correction setting, a
    %   logical (true/false), which enabled or disables the use of
    %   specialized LR_* functions that have additional height-thresholding
    %   and higher-order correction factors built into the data correction.
    
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
        fitLog logical
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
        
        function obj = ViscoFit(forces,times,indentations,tipSize,varargin)
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
                obj.tipGeom = "spherical";
                obj.fitLog = false;
                obj.minTimescale = 1e-4; % Default
                if ~isempty(varargin)
                    if isa(varargin{1},'struct')
                        % Grab the settings structure
                        inputSettings = varargin{1};
                        
                        % Store the settings in the class
                        if ~logical(inputSettings.fitLog)
                            temp = cellfun(@(nu,t) nu.*ones(size(t)),inputSettings.nu,times,'UniformOutput',false);
                            obj.nu = horzcat(temp{:});
                            obj.nu_cell = temp;
                            
                            if logical(inputSettings.thinSample)
%                                 temp = cellfun(@(h,t) h.*ones(size(t)),inputSettings.pixelHeight,times,'UniformOutput',false);
%                                 obj.pixelHeight = horzcat(temp{:});
                                obj.pixelHeight_cell = inputSettings.pixelHeight;
                            end
                        else
                            tempdt = cellfun(@(t) round(mode(gradient(t)),1,'significant'),times,'UniformOutput',false);
                            temp = cellfun(@(x,t,dt)log_scale(x.*ones(size(t)),t,mode(dt),t(end)),inputSettings.nu,times,tempdt,'UniformOutput',false);
                            obj.nu = horzcat(temp{:});
                            obj.nu_cell = temp;
                            
                            if logical(inputSettings.thinSample)
%                                 tempdt = cellfun(@(t) round(mode(gradient(t)),1,'significant'),times,'UniformOutput',false);
%                                 temp = cellfun(@(x,t,dt)log_scale(x.*ones(size(t)),t,mode(dt),t(end)),inputSettings.pixelHeight,times,tempdt,'UniformOutput',false);
%                                 obj.pixelHeight = horzcat(temp{:});
                                obj.pixelHeight_cell = inputSettings.pixelHeight;
                            end
                        end
                        obj.tipGeom = string(inputSettings.tipGeom);
                        obj.fitLog = inputSettings.fitLog;
                        obj.minTimescale = inputSettings.minTimescale;
                        
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
                        error('You are not passing the settings correctly to the ViscoFit Class Initialization. Please ensure the fifth argument is a structure containing your settings.');
                    end
                    
                end
                
                if ~obj.fitLog
                    % Use the full-fidelity data (no log sampling)
                    obj.forces = horzcat(forces{:});
                    obj.forces_cell = forces;

                    obj.times = horzcat(times{:});
                    obj.times_cell = times;

                    obj.indentations = horzcat(indentations{:});
                    obj.indentations_cell = indentations;

                    temp = cellfun(@(t) round(mode(gradient(t)),1,'significant').*ones(size(t)),times,'UniformOutput',false);
                    obj.dts = horzcat(temp{:});
                    obj.dts_cell = temp;

                    temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,times,'UniformOutput',false);
                    obj.tipSize = horzcat(temp{:});
                    obj.tipSize_cell = temp;
                else
                    % Log-sample all of the data before performing any
                    % operations. Note: this reduces the data quality, but
                    % drastically improves the fitting speed for legacy
                    % methods (NLS, in particular).
                    tempdt = cellfun(@(t) round(mode(gradient(t)),1,'significant'),times,'UniformOutput',false);
                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),times,times,tempdt,'UniformOutput',false);
                    obj.times = horzcat(temp{:});
                    obj.times_cell = temp;
                    
                    temp = cellfun(@(t,x) round(mode(gradient(t)),1,'significant').*ones(size(x)),times,obj.times_cell,'UniformOutput',false);
                    obj.dts = horzcat(temp{:});
                    obj.dts_cell = temp;
                    
                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),forces,times,obj.dts_cell,'UniformOutput',false);
                    obj.forces = horzcat(temp{:});
                    obj.forces_cell = temp;

                    temp = cellfun(@(x,t,dt)log_scale(x,t,mode(dt),t(end)),indentations,times,obj.dts_cell,'UniformOutput',false);
                    obj.indentations = horzcat(temp{:});
                    obj.indentations_cell = temp;
                    
                    temp = cellfun(@(r,t) r.*ones(size(t)),tipSize,obj.times_cell,'UniformOutput',false);
                    obj.tipSize = horzcat(temp{:});
                    obj.tipSize_cell = temp;
                    
                end
                
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
        
        function fitStruct = fitData(obj, varargin)
            % FITDATA Fit a Viscoelastic Model to the Class Data 
            %   This function takes in a variety of settings in addition to
            %   the data already provided to the class and performs an
            %   optimization procedure based on those settings. Defaults
            %   are set below, and the order in which they are provided
            %   must be adhered to (unless the user manually reprograms
            %   this class). In particular, the number of iterations
            %   (n_iterations) may need to be increased, as with the number
            %   of solver iterations (n_fitIterations) and maximum number
            %   of elements in the viscoelastic series (n_elements). This
            %   function will iteratively introduce viscoelastic elements,
            %   feeding forward the results from the previous iteration,
            %   until the optimization has been performed for a model
            %   configuration with a number of viscoelastic elements 
            %   equal to n_elements.
            
            % Initialize Output Structure
            fitStruct = struct;
            
            % Default Settings
            solver = 'nelder-mead';     % Fit using Nelder-Mead Simplex
            model = 'maxwell';          % Use Generalized Maxwell Model
            n_elements = 3;             % Fit iteratively for up to 3 elements
            elasticSetting = 1;         % Include Elastic Term
            fluidSetting = 0;           % No Steady-State Fluidity
            n_iterations = 100;         % Use 100 random initializations as a default
            n_fitIterations = 1e4;      % No. of iterations for solver
            errortype = 'sse';          % Error model to use
            N_workers = [];             % Number of workers for parpool

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
                    errortype = fitOpts.errortype;
                end
                try
                    N_workers = fitOpts.N_workers;
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
            fitStruct.errortype = errortype;
            fitStruct.ViscoClass = obj;
            
            % Get the correct objective function for optimization
            switch model
                case 'maxwell'
                    objFunc = @obj.SSE_Maxwell;
                case 'voigt'
                    objFunc = @obj.SSE_Voigt;
                case 'plr'
                    objFunc = @obj.SSE_PLR;
                    % The storage roster for plr is different, and requires
                    % the second index to maintain consistency. Thus, we
                    % have to manually force the second term to be fit by
                    % including the fluidity. This second position actually
                    % corresponds to the exponent, alpha.
                    fluidSetting = 1;
                    fitStruct.fluidSetting = fluidSetting;
                case 'custom'
                    objFunc = @obj.customFunc;
                otherwise
                    error('Your chosen solver-model combination is not implemented yet.');
            end
            
            % Create placeholders for the data we will obtain from our
            % optimization attempts
            fitStruct.bestParams = {};
            fitStruct.paramPopulation = {};
            fitStruct.paramPopulationResiduals = {};
            fitStruct.elasticFitTime = {};
            fitStruct.fitTime = {};
            fitStruct.upperParamCI = {};
            fitStruct.lowerParamCI = {};
            
            if ~isempty(gcp('nocreate'))
               % Get the current pool
                poolobj = gcp('nocreate');
                delete(poolobj);
            end
            
            % Make a fresh pool
            if isempty(N_workers)
                poolobj = parpool('IdleTimeout', 120);
            else
                poolobj = parpool(N_workers,'IdleTimeout', 120);
            end

            % Send the class to the workers
            addAttachedFiles(poolobj, {'ViscoFit.m','LR_Maxwell.m','LR_Voigt.m','LR_PLR.m'})

            % Start the timer
            tic;
            
            % Begin the iterative term introduction loop
            for i = 1:n_elements   
                
                % Look to see if there are old results available to provide
                % intelligent guesses for our initial conditions. This will
                % only occur for iterations beyond the first
                if i > 1
                    % Our current loop has an array of parameters
                    % two-larger than the previous (for the generalized
                    % spring-dashpot rheology models)
                    beta_in = NaN(length(fitStruct.bestParams{i-1})+2,1);
                    
                    % Overwrite the NaN values with the previous optimal
                    % parameters.
                    beta_in(1:length(fitStruct.bestParams{i-1})) = fitStruct.bestParams{i-1};
                else
                    % For the first iteration, there will be four
                    % parameters in total: the elastic element (1), the
                    % steady-state fluidity (2), and two parameters for the
                    % first viscoelastic element (3,4). Note, based on the
                    % elasticSetting and fluidSetting provided to this
                    % function, the values in (1) and (2) may never be
                    % updated
                    if strcmpi(model,'plr')
                        beta_in = NaN(2,1);
                    else
                        beta_in = NaN(4,1);
                    end
                end
                
                % Make our storage variables
                beta_dist = zeros(length(beta_in),n_iterations);
                beta_dist_elastic = zeros(1,n_iterations);
                beta0_dist = zeros(size(beta_dist));
                residual_dist = NaN(1,n_iterations);
                residual_dist_elastic = NaN(size(residual_dist));
                
                % Make the upper and lower bounds for this model
                [tauInds,modulusInds] = getParamIndices(beta_in);
                ub = zeros(size(beta_in))+eps;
                ub_rand = [];
                lb = zeros(size(beta_in));
                lb_rand = [];
                
                switch model
                    case 'maxwell'
                        % Moduli are limited to "reasonable" bounds for
                        % viscoelastic materials
                        ub(modulusInds) = 1e12;
                        lb(modulusInds) = 1e-2;

                        tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                        ub(tauInds) = tauCenters*10;
                        lb(tauInds) = tauCenters/10;
                        
                        if ~elasticSetting
                            ub(1) = eps;
                            lb(1) = 0;
                        end

                        if fluidSetting
                            ub(2) = max(tauCenters)*1e2;
                            lb(2) = 10^( floor(min(log10(obj.dts)))+1 );
                        end
                        
                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand(modulusInds) = 6;
                        ub_rand(tauInds) = log10(ub(tauInds));
                        lb_rand(modulusInds) = 1;
                        lb_rand(tauInds) = log10(lb(tauInds));
                        
                    case 'voigt'
                        % Compliances are limited to "reasonable" bounds
                        % for viscoelastic materials
                        ub(modulusInds) = 1e2;
                        lb(modulusInds) = 1e-12;

                        tauCenters = obj.minTimescale.*(10.^( (1:length(ub(3:2:end)))-1 ));
                        ub(tauInds) = tauCenters*10;
                        lb(tauInds) = tauCenters/10;
                        
                        if ~elasticSetting
                            ub(1) = eps;
                            lb(1) = 0;
                        end

                        if fluidSetting
                            ub(2) = 1;
                            lb(2) = 0;
                        end

                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand(modulusInds) = -1;
                        ub_rand(tauInds) = log10(ub(tauInds));
                        lb_rand(modulusInds) = -6;
                        lb_rand(tauInds) = log10(lb(tauInds));
                        
                    case 'plr'
                        % Power Law Rheology Roster:
                        % [E_0 alpha]
                        ub = [1e12 1];
                        lb = [1e-2 0];
                        
                        % Restrict the range of random guesses, if desired.
                        % Otherwise, they should be set equal to ub & lb
                        ub_rand = log10(ub);
                        lb_rand = [log10(lb(1)) -3];
                        
                    case 'custom'
                        
                        % Define the upper and lower bounds for your custom
                        % function here. The output from this region should
                        % be two arrays, ub and lb, which are the same
                        % length as the number of parameters in the model.
                        % ub = [...];
                        % lb = [...];
                        
                end
                
                preFitting = toc;
                switch solver
                    case 'nelder-mead'
                        
                        options = optimset('Display','none',...
                                    'PlotFcns',[],...
                                    'MaxFunEvals',n_fitIterations,...
                                    'MaxIter',n_fitIterations,...
                                    'TolFun',1e-60,...
                                    'TolX',1e-60);
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,options);
                            end
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,options);
                        end
                            
                    case 'annealing'
                        
                        nelderopts = optimset('Display','none',...
                                        'PlotFcns',[],...
                                        'MaxFunEvals',n_fitIterations,...
                                        'MaxIter',n_fitIterations,...
                                        'TolFun',1e-60,...
                                        'TolX',1e-60);
                                
                        annealopts = struct(...
                                        'CoolSched',@(T) (.4*T),...
                                        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                        'InitTemp',1,...
                                        'MaxConsRej',1000,...
                                        'MaxSuccess',20,...
                                        'MaxTries',300,...
                                        'StopTemp',1e-5,...
                                        'StopVal',-Inf,...
                                        'Verbosity',0);
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fminsearch(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,nelderopts);
                            end
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = annealOpt(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,annealopts,nelderopts);
                        end
                        
                    case 'nls'
                        
                        fminoptions = optimoptions('fmincon','Algorithm','sqp',...
                                        'MaxFunctionEvaluations', n_fitIterations,...
                                        'MaxIterations', n_fitIterations,...
                                        'FiniteDifferenceType','central',...
                                        'FunctionTolerance', 1e-60,...
                                        'OptimalityTolerance', 1e-60,...
                                        'Display', 'none');
                        
                        if i == 1 && elasticSetting
                            % Fit the elastic term separately for the first
                            % iteration. Future iterations have the "best
                            % fit" elastic term included from the prior
                            % optimization attempt
                            
                            % Clock the timer
                            preElasticTime = toc;
                            
                            parfor j = 1:n_iterations
                                % Get the grid search starting position
                                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{j});
                                [beta_dist_elastic(j),residual_dist_elastic(j)] = fmincon(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
                            end
                            
                            % Clock the timer and save the fitting time
                            postElasticTime = toc;
                            fitStruct.elasticFitTime = horzcat(fitStruct.elasticFitTime,{postElasticTime-preElasticTime});
                            
                            % Find the best elastic parameter
                            [~,idx] = min(residual_dist_elastic,[],'omitnan');
                            beta_in(1) = beta_dist_elastic(:,idx);
                        end
                        
                        % See which parameters are new this time, so that
                        % information can be fed to our
                        % random-guess-generation function
                        newInds = isnan(beta_in);
                        
                        parfor j = 1:n_iterations
                            % Get the grid search starting position
                            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
                            beta0_dist(:,j) = beta0;
                            [beta_dist(:,j),residual_dist(j)] = fmincon(@(x)objFunc(x,elasticSetting,fluidSetting,errortype),beta0,[],[],[],[],lb,ub,[],fminoptions);
                        end
                        
                    otherwise
                        error('That solver is not supported.')
                end
                postFitting = toc;
                
                % Find the best-fit parameters from our population
                [~,idx] = min(residual_dist,[],'omitnan');
                if size(idx,2)>1 idx = (idx(1)); end
                
                % Store the best-fit parameters to be fed forward, and also
                % save the "less optimal" sets for statistical treatment
                % later
                fprintf('Optimal Parameters, Iteration %d:\n',i)
                disp(beta_dist(:,idx));
                fitStruct.bestParams = horzcat(fitStruct.bestParams,{beta_dist(:,idx)});
                fitStruct.paramPopulation = horzcat(fitStruct.paramPopulation,{beta_dist});
                fitStruct.paramPopulationResiduals = horzcat(fitStruct.paramPopulationResiduals,{residual_dist});
                [tempub,templb] = obj.getParamsCI(beta_dist,0.95);
                fitStruct.upperParamCI = horzcat(fitStruct.upperParamCI,{tempub});
                fitStruct.lowerParamCI = horzcat(fitStruct.lowerParamCI,{templb});
                
                % Store the timing for this model configuration fit
                fitStruct.fitTime = horzcat(fitStruct.fitTime,{postFitting-preFitting});
                
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
                
            end % End Iterative Term Introduction Loop
                        
        end % End fitData()
        
    end % End Methods
    
    methods (Static = true)
        
        % Empty
        
    end % End Static Methods
    
end

