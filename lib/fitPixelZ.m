function [relaxanceMap_out,retardanceMap_out,alphaMap_out,frequencyMap_out,bestParamsMap_out,paramPopulationMap_out,paramPopulationResidualsMap_out,elasticFitTimeMap_out,fitTimeMap_out] = fitPixelZ(n_terms,n_iterations,model,solver,n_fitIterations,minTimescaleTemp,dataInDistrib,paramLen,objFuncMap,elasticSetting,fluidSetting,bestParamsMap,smoothOpt,windowsize,thinPixel)
%FITPIXELZ Perform a Viscoelastic Fit to One Pixel in a QI Map using the
%Z-Transform Method
%   This function takes in a variety of settings in addition to
%   the data already provided to the class and performs an
%   optimization procedure based on those settings.
%
%   This function will result in an optimized parameter set
%   for a single pixel in a force map, and is called by the "fitMap..."
%   function family (depending on which approach is being used).

% Look to see if there are old results available to provide
% intelligent guesses for our initial conditions. This will
% only occur for iterations beyond the first
if n_terms > 1
    % Our current loop has an array of parameters
    % two-larger than the previous (for the generalized
    % spring-dashpot rheology models)
    beta_in = NaN(paramLen+2,1);

    % Overwrite the NaN values with the previous optimal
    % parameters.
    beta_in(1:paramLen) = bestParamsMap;
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

% Convert datatypes
if isdistributed(dataInDistrib)
    dataIn = gather(dataInDistrib);
elseif numel(dataInDistrib) == 1
    dataIn = [dataInDistrib{:}];
end

% Prevent having an output that is not set
elasticFitTimeMap_out = {0};

switch model
    case 'maxwell'
        % Moduli are limited to "reasonable" bounds for
        % viscoelastic materials
        ub(modulusInds) = 1e12;
        lb(modulusInds) = 1e-2;

        % Use minTimescale
        tauCenters = minTimescaleTemp.*(10.^( (1:length(ub(3:2:end)))-1 ));
%         ub(tauInds) = tauCenters*10;
%         lb(tauInds) = tauCenters/10;
        
        % Don't use minTimescale
        lb(tauInds) = 10^(ceil(log10(mode(dataIn{2})))-2);
        ub(tauInds) = 1./lb(tauInds);

        if ~elasticSetting
            ub(1) = eps;
            lb(1) = 0;
        end

        if fluidSetting
            ub(2) = max(tauCenters)*1e2;
            lb(2) = lbFluidity;
        end

        % Restrict the range of random guesses, if desired.
        % Otherwise, they should be set equal to ub & lb
        ub_rand(modulusInds) = 6;
        ub_rand(tauInds) = ceil(log10(max(dataIn{1})))+1;
        lb_rand(modulusInds) = 0;
        lb_rand(tauInds) = floor(log10(dataIn{8}))-1;

    case 'voigt'
        % Compliances are limited to "reasonable" bounds
        % for viscoelastic materials
        ub(modulusInds) = 1e2;
        lb(modulusInds) = 1e-12;

        tauCenters = minTimescaleTemp.*(10.^( (1:length(ub(3:2:end)))-1 ));
%         ub(tauInds) = tauCenters*10;
%         lb(tauInds) = tauCenters/10;

        % Don't use minTimescale
        lb(tauInds) = 10^(ceil(log10(mode(dataIn{2})))-2);
        ub(tauInds) = 1./lb(tauInds);

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
        ub_rand(modulusInds) = 0;
        ub_rand(tauInds) = ceil(log10(max(dataIn{1})))+1;
        lb_rand(modulusInds) = -6;
        lb_rand(tauInds) = floor(log10(dataIn{8}))-1;

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

if any(isnan(dataIn{3}))
    % This is the "bad pixel" trigger! Skip this pixel
    % and enter in NaN for the output data.
    elasticFitTimeMap_out = {0};
    fitTimeMap_out = {0};

    relaxanceMap_out = NaN;
    retardanceMap_out = NaN;
    alphaMap_out = NaN;
    frequencyMap_out = NaN;
    bestParamsMap_out = {beta_dist(:,1)};
    paramPopulationMap_out = {beta_dist};
    paramPopulationResidualsMap_out = {residual_dist};
    return;
end

% Get Z-Transform Data from Curve
% Do NOT save smoothed data. Only use the smoothing for fitting!
[Q_hz,F_t,~,h_t,~,f_hz,alphaInit] = zTransformCurve(dataIn,'none',windowsize,thinPixel);

U_hz = 1./(Q_hz);
relaxanceMap = Q_hz;
retardanceMap = U_hz;
alphaMap = alphaInit;
frequencyMap = f_hz;

tic;
switch solver
    case 'nelder-mead'

        options = optimset('Display','iter-detailed',...
                    'PlotFcns',[],...
                    'MaxFunEvals',n_fitIterations,...
                    'MaxIter',n_fitIterations,...
                    'TolFun',0,...
                    'TolX',0);

        if n_terms == 1 && elasticSetting
            % Fit the elastic term separately for the first
            % iteration. Future iterations have the "best
            % fit" elastic term included from the prior
            % optimization attempt

            % Clock the timer
            preElasticTime = toc;

            for k = 1:n_iterations
                % Get the grid search starting position
                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                [beta_dist_elastic(k),residual_dist_elastic(k)] = fminsearch(@(x)objFuncMap(dataIn,x,ub(1),lb(1),smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,options);
            end

            % Clock the timer and save the fitting time
            postElasticTime = toc;
            elasticFitTimeMap_out = {postElasticTime-preElasticTime};

            % Find the best elastic parameter
            [~,idx] = min(residual_dist_elastic,[],'omitnan');
            beta_in(1) = beta_dist_elastic(:,idx);
            
            % Adjust random starting range
            delt = abs(ub_rand(1)-lb_rand(1))/2;
            ub_rand(modulusInds) = round(log10(beta_in(1)))+delt;
            lb_rand(modulusInds) = round(log10(beta_in(1)))-delt;
            
        end

        % See which parameters are new this time, so that
        % information can be fed to our
        % random-guess-generation function
        newInds = isnan(beta_in);
        
        for k = 1:n_iterations
            % Get the grid search starting position
            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
            beta0_dist(:,k) = beta0;
            [beta_dist(:,k),residual_dist(k)] = fminsearch(@(x)objFuncMap(dataIn,x,ub,lb,smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,options);
        end
        
%         % Observe Smoothing Results!
%         try
%             figure(hfig)
%             clf;
%         catch
%             hfig = figure;
%         end
%         tiledlayout(1,4);
%         
%         % Find the best-fit parameters from our population
%         [~,idx] = min(residual_dist,[],'omitnan');
%         if size(idx,2)>1
%             idx = (idx(1));
%         end
%         
%         omega = linspace(-pi,pi,numel(Q_hz));
%         omega = omega(omega>=0);
%         [Q_func,Q_storage,Q_loss] = getZModels('maxwell',n_terms,alphaInit);
%         
%         plotParams = beta_dist(:,idx);
% %         plotParams = [1e4 0 1e4 5e-4];
%         
%         nexttile
%         plot(dataIn{4},dataIn{3},'b-','linewidth',3)
%         hold on
%         plot(h_t,F_t,'r-','linewidth',3)
% %         plot(h_t,LR_Maxwell(beta0_dist(:,idx),dataIn{1},dataIn{2},...
% %             dataIn{4},dataIn{5},dataIn{6},dataIn{7},elasticSetting,...
% %             fluidSetting,dataIn{9},dataIn{10}),'g-','linewidth',2);
%         plot(h_t,LR_Maxwell(plotParams,dataIn{1},dataIn{2},dataIn{4},...
%             dataIn{5},dataIn{6},dataIn{7},elasticSetting,fluidSetting,...
%             dataIn{9},dataIn{10}),'k-','linewidth',2);
%         xlabel('Indentation [m]')
%         ylabel('Force [N]')
% %         legend('Original','Corrected','Initial Fit','Final Fit','location','best')
%         legend('Original','Corrected','Final Fit','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(f_hz(f_hz>=0),Q_hz(f_hz>=0),'bo','linewidth',3)
%         hold on
%         plot(f_hz((end-numel(omega)+1):end),Q_func(beta0_dist(:,idx),omega),'b-','linewidth',3)
%         plot(f_hz((end-numel(omega)+1):end),Q_func(plotParams,omega),'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         legend('Data','Initial Fit','Final Fit','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(f_hz(f_hz>=0),real(Q_hz(f_hz>=0)),'bo','linewidth',3)
%         hold on
%         plot(f_hz((end-numel(omega)+1):end),Q_storage(beta0_dist(:,idx),omega),'b-','linewidth',3)
%         plot(f_hz((end-numel(omega)+1):end),Q_storage(plotParams,omega),'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Storage Modulus [Pa]')
%         legend('Data','Initial Fit','Final Fit','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(f_hz(f_hz>=0),abs(imag(Q_hz(f_hz>=0))),'bo','linewidth',3)
%         hold on
%         plot(f_hz((end-numel(omega)+1):end),Q_loss(beta0_dist(:,idx),omega),'b-','linewidth',3)
%         plot(f_hz((end-numel(omega)+1):end),Q_loss(plotParams,omega),'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Loss Modulus [Pa]')
%         legend('Data','Initial Fit','Final Fit','location','best')
%         grid on
%         hold off
%         
%         disp('pause');
        
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

        if n_terms == 1 && elasticSetting
            % Fit the elastic term separately for the first
            % iteration. Future iterations have the "best
            % fit" elastic term included from the prior
            % optimization attempt

            % Clock the timer
            preElasticTime = toc;

            for k = 1:n_iterations
                % Get the grid search starting position
                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                [beta_dist_elastic(k),residual_dist_elastic(k)] = fminsearch(@(x)objFuncMap(dataIn,x,ub(1),lb(1),smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,nelderopts);
            end

            % Clock the timer and save the fitting time
            postElasticTime = toc;
            elasticFitTimeMap_out = {postElasticTime-preElasticTime};

            % Find the best elastic parameter
            [~,idx] = min(residual_dist_elastic,[],'omitnan');
            beta_in(1) = beta_dist_elastic(:,idx);
        end

        % See which parameters are new this time, so that
        % information can be fed to our
        % random-guess-generation function
        newInds = isnan(beta_in);

        for k = 1:n_iterations
            % Get the grid search starting position
            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
            beta0_dist(:,k) = beta0;
            [beta_dist(:,k),residual_dist(k)] = annealOpt(@(x)objFuncMap(dataIn,x,ub,lb,smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,annealopts,nelderopts);
        end

    case 'nls'

        fminoptions = optimoptions('fmincon','Algorithm','sqp',...
                        'MaxFunctionEvaluations', n_fitIterations,...
                        'MaxIterations', n_fitIterations,...
                        'FiniteDifferenceType','central',...
                        'FunctionTolerance', 0,...
                        'OptimalityTolerance', 1e-30,...
                        'StepTolerance', 1e-10,...
                        'Display', 'iter-detailed');

        if n_terms == 1 && elasticSetting
            % Fit the elastic term separately for the first
            % iteration. Future iterations have the "best
            % fit" elastic term included from the prior
            % optimization attempt

            % Clock the timer
            preElasticTime = toc;

            for k = 1:n_iterations
                % Get the grid search starting position
                beta0 = getfield(logspace(ub_rand(1),lb_rand(1),n_iterations),{k});
                [beta_dist_elastic(k),residual_dist_elastic(k)] = fmincon(@(x)objFuncMap(dataIn,x,ub(1),lb(1),smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,[],[],[],[],lb(1),ub(1),[],fminoptions);
            end

            % Clock the timer and save the fitting time
            postElasticTime = toc;
            elasticFitTimeMap_out = {postElasticTime-preElasticTime};

            % Find the best elastic parameter
            [~,idx] = min(residual_dist_elastic,[],'omitnan');
            beta_in(1) = beta_dist_elastic(:,idx);
        end

        % See which parameters are new this time, so that
        % information can be fed to our
        % random-guess-generation function
        newInds = isnan(beta_in);

        for k = 1:n_iterations
            % Get the grid search starting position
            beta0 = makeRandomParams(beta_in,ub_rand,lb_rand,elasticSetting,fluidSetting,newInds);
            beta0_dist(:,k) = beta0;
            [beta_dist(:,k),residual_dist(k)] = fmincon(@(x)objFuncMap(dataIn,x,ub,lb,smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting),beta0,[],[],[],[],lb,ub,[],fminoptions);
        end
        
    otherwise
        error('That solver is not supported.')
        
end
postFitting = toc;

% Find the best-fit parameters from our population
[~,idx] = min(residual_dist,[],'omitnan');
if size(idx,2)>1
    idx = (idx(1));
end

% Store the best-fit parameters to be fed forward, and also
% save the "less optimal" sets for statistical treatment
% later
relaxanceMap_out = {relaxanceMap};
retardanceMap_out = {retardanceMap};
alphaMap_out = {alphaMap};
frequencyMap_out = {frequencyMap};
bestParamsMap_out = {beta_dist(:,idx)};
paramPopulationMap_out = {beta_dist};
paramPopulationResidualsMap_out = {residual_dist};

% Store the timing for this model configuration fit
fitTimeMap_out = {postFitting};

% Release Variables Manually
dataIn = [];
beta_dist = [];
beta_dist_elastic = [];
beta0_dist = [];
residual_dist = [];
residual_dist_elastic = [];
tauInds = [];
modulusInds = [];
ub = [];
lb = [];

end

