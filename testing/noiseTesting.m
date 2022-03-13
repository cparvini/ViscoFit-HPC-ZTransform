clear all
close all
clc

%% Simulated Data Noise Testing

% Select a simulation file
simFileFolder = uigetdir(pwd,'Select a Directory Containing a Dataset .mat File');

% Check to see if there are subdirectories
dirContents = dir(simFileFolder);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root filesep sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
else
    Folders = {simFileFolder};
end

outStruct = struct();

for i_dir = 1:length(Folders)
    
    % Start with a clean slate
    clearvars -except i_dir Folders simFileFolder outStruct
    close all
    clc
    fprintf('Analyzing Directory #%d of %d\n',i_dir,length(Folders));

    path = Folders{i_dir};
    
    % Load Sim Data
    settingsFile = dir([path filesep '*.mat']);
    settingsFile(~endsWith({settingsFile.name},'Settings.mat')) = [];

    dataFile = dir([path filesep '*.mat']);
    dataFile(endsWith({dataFile.name},'Settings.mat')) = [];

    % Settings
    n = 6;
    SNR = 50;
    n_terms = 1;
    n_iterations = 100;
    model = 'maxwell';
    solver = 'nelder-mead';
    n_fitIterations = 3e3;
    paramLen = 2*n_terms + 2;
    smoothOpt = 'none';
    windowsize = 0.05;
    thinPixel = false;

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

    % Load variables from the file
    simStruct = load(fullfile(dataFile.folder,dataFile.name));
    load(fullfile(settingsFile.folder,settingsFile.name));

    F = simStruct.F;
    z = simStruct.z;
    d = simStruct.d;
    h = -(z-d);
    t = simStruct.t;
    dt = settingsStruct.dt;
    r_tip = settingsStruct.r_tip;
    nu = settingsStruct.nu_sample;
    elasticSetting = (settingsStruct.E0 > 0);
    fluidSetting = (settingsStruct.phi_f > 0);
    trueParams = [settingsStruct.E0;settingsStruct.phi_f];
    for j = 1:size(settingsStruct.linearParamData,1)
        trueParams = vertcat(trueParams,...
            [settingsStruct.linearParamData(j,2);...
            settingsStruct.linearParamData(j,4)]);
    end

    clearvars simStruct settingsStruct

    % Trim data
    [~,zid] = max(-z);
    hid = find(h>=0,1);
    F = F(hid:zid);
    h = h(hid:zid);
    t = t(hid:zid);
    t = t - t(1);

    % Remove the 0-time point if it exists
    F(t==0) = [];
    h(t==0) = [];
    d(t==0) = [];
    z(t==0) = [];
    t(t==0) = [];

    if ~isrow(F)
        F = F';
    end
    if ~isrow(h)
        h = h';
    end
    if ~isrow(t)
        t = t';
    end
    if ~isrow(z)
        z = z';
    end
    if ~isrow(d)
        d = d';
    end
    fs = 1/dt;
    
    % Create the "True" dataIn
    dataInTrue = cell(1,10);
    dataInTrue{1} = t;                   % time (s)
    dataInTrue{2} = dt*ones(size(t));    % timestep (s)
    dataInTrue{3} = F;                   % force (N)
    dataInTrue{4} = h;                   % indentation (m)
    dataInTrue{5} = r_tip*ones(size(t)); % tip radius (m)
    dataInTrue{6} = nu*ones(size(t));    % sample Poisson's Ratio (unitless)
    dataInTrue{7} = "spherical";         % tip geometry
    dataInTrue{8} = 10^floor(log10(dt)); % minimum timescale for fitting (s)
    dataInTrue{9} = thinPixel;           % (TS) thin sample correction T/F
    dataInTrue{10} = NaN;                % sample finite thickness (for TS)

    % Add the desired noise
    if n > 6
        error('You have asked for too many noise types. The max is six (6).');
    end
    dataIn = cell(n,1);

    for j = 1:n
        dataIn{j} = cell(1,10);
        dataIn{j}{1} = t;                   % time (s)
        dataIn{j}{2} = dt*ones(size(t));    % timestep (s)
        dataIn{j}{4} = h;                   % indentation (m)
        dataIn{j}{5} = r_tip*ones(size(t)); % tip radius (m)
        dataIn{j}{6} = nu*ones(size(t));    % sample Poisson's Ratio (unitless)
        dataIn{j}{7} = "spherical";         % tip geometry
        dataIn{j}{8} = 10^floor(log10(dt)); % minimum timescale for fitting (s)
        dataIn{j}{9} = thinPixel;           % (TS) thin sample correction T/F
        dataIn{j}{10} = NaN;                % sample finite thickness (for TS)

        signalPower = sum(F.^2)/numel(F);

        switch j
            case 1
                % white gaussian noise
                dataIn{j}{3} = awgn(F,SNR,'measured');

            case 2
                temp = dsp.ColoredNoise('pink',numel(F),1);
                noise = temp();
                if ~isrow(noise)
                    noise = noise';
                end
                noisePower = sum(noise.^2)/numel(noise);
                scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
                dataIn{j}{3} = F + noise*scaleFactor;

            case 3
                temp = dsp.ColoredNoise('white',numel(F),1);
                noise = temp();
                if ~isrow(noise)
                    noise = noise';
                end
                noisePower = sum(noise.^2)/numel(noise);
                scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
                dataIn{j}{3} = F + noise*scaleFactor;

            case 4
                temp = dsp.ColoredNoise('brown',numel(F),1);
                noise = temp();
                if ~isrow(noise)
                    noise = noise';
                end
                noisePower = sum(noise.^2)/numel(noise);
                scaleFactor = sqrt(signalPower./(noisePower*(10^(SNR/10))));
                dataIn{j}{3} = F + noise*scaleFactor;

            case 5
                temp = dsp.ColoredNoise('blue',numel(F),1);
                noise = temp();
                if ~isrow(noise)
                    noise = noise';
                end
                noisePower = sum(noise.^2)/numel(noise);
                scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
                dataIn{j}{3} = F + noise*scaleFactor;

            case 6
                temp = dsp.ColoredNoise('purple',numel(F),1);
                noise = temp();
                if ~isrow(noise)
                    noise = noise';
                end
                noisePower = sum(noise.^2)/numel(noise);
                scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
                dataIn{j}{3} = F + noise*scaleFactor;

            otherwise
                continue;
        end

    %     % Quick Look
    %     plot(t,F)
    %     hold on
    %     plot(t,dataIn{j}{3})
    %     legend('Original','Noisy')
    %     hold off

    end

    % Pre-allocation
    relaxanceMap = cell(n,1);
    retardanceMap = cell(n,1);
    alphaMap = cell(n,1);
    frequencyMap = cell(n,1);
    bestParamsMap = cell(n,1);
    paramPopulationMap = cell(n,1);
    paramPopulationResidualsMap = cell(n,1);
    elasticFitTimeMap = cell(n,1);
    fitTimeMap = cell(n,1);

    parfor j = 1:n
        bestParamsMap{j} = NaN(4,1);
        paramPopulationMap{j} = NaN(4,n_iterations);
        paramPopulationResidualsMap{j} = NaN(1,n_iterations);
        elasticFitTimeMap{j} = NaN;
        fitTimeMap{j} = NaN;
        relaxanceMap{j} = NaN(n,1);
        retardanceMap{j} = NaN(n,1);
        alphaMap{j} = NaN;
        frequencyMap{j} = NaN(n,1);
    end
    
    % Kill parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);

    % Perform the fitting
    for j = 1:n

        fprintf('Processing Iteration %d\n',j);

        [relaxanceMap(j),...
            retardanceMap(j),...
            alphaMap(j),...
            frequencyMap(j),...
            bestParamsMap(j),...
            paramPopulationMap(j),...
            paramPopulationResidualsMap(j),...
            elasticFitTimeMap(j),...
            fitTimeMap(j)] = fitPixelZ(n_terms,n_iterations,model,solver,n_fitIterations,dataIn{j}{8},dataIn(j),paramLen,objFuncMap,elasticSetting,fluidSetting,bestParamsMap(j),smoothOpt,windowsize,thinPixel);

        fprintf('Total Time: %4.1f min; Elastic Time: %4.1f min\n\n',fitTimeMap{j}/60,elasticFitTimeMap{j}/60);

    end
    
    % Get Z-Transform Data from Noise-free Curve
    [Q_hz_true,~,~,~,~,f_hz_true,~] = zTransformCurve(dataInTrue,'none',windowsize,thinPixel);

    for j = 1:n

        % Plot fit results
        try
            figure(hfig)
            clf;
        catch
            hfig = figure;
            hfig.Position = [0 25 1600 600];
        end
        tiledlayout(1,4,'TileSpacing','compact');

        % Labels
        switch j
            case 1
                plotTitle = ['White Gaussian Noise, SNR ' num2str(SNR)];
                saveLabel = 'WhiteGaussian';
            case 2
                plotTitle = ['Pink Noise, SNR ' num2str(SNR)];
                saveLabel = 'PinkNoise';
            case 3
                plotTitle = ['White Noise, SNR ' num2str(SNR)];
                saveLabel = 'WhiteNoise';
            case 4
                plotTitle = ['Brown Noise, SNR ' num2str(SNR)];
                saveLabel = 'BrownNoise';
            case 5
                plotTitle = ['Blue Noise, SNR ' num2str(SNR)];
                saveLabel = 'BlueNoise';
            case 6
                plotTitle = ['Purple Noise, SNR ' num2str(SNR)];
                saveLabel = 'PurpleNoise';
            otherwise
                continue;
        end

        % Grab the data
        plotParams = bestParamsMap{j};
        Q_hz = relaxanceMap{j};
        f_hz = frequencyMap{j};
        alphaInit = alphaMap{j};

        omega = linspace(-pi,pi,numel(Q_hz));
        omega = omega(omega>=0);
        [Q_func,Q_storage,Q_loss] = getZModels('maxwell',n_terms,alphaInit);
        
        nexttile
        plot(dataIn{j}{4},dataIn{j}{3},'b-','linewidth',7)
        hold on
        plot(h,F,'r-','linewidth',3)
        xlabel('Indentation [m]')
        ylabel('Force [N]')
        legend('Fitting Data','Truth','location','best')
        grid on
        title(plotTitle)
        hold off

        nexttile
        scatter(f_hz(f_hz>=0),Q_hz(f_hz>=0),'bo','linewidth',3)
        hold on
        plot(f_hz((end-numel(omega)+1):end),Q_func(plotParams,omega),'r-','linewidth',7)
        plot(f_hz_true((end-numel(omega)+1):end),Q_hz_true((end-numel(omega)+1):end),'g-','linewidth',3)
        xlabel('Frequency [Hz]')
        ylabel('Relaxance [Pa]')
        legend('Data','Fit','Truth','location','best')
        grid on
        hold off

        nexttile
        scatter(f_hz(f_hz>=0),real(Q_hz(f_hz>=0)),'bo','linewidth',3)
        hold on
        plot(f_hz((end-numel(omega)+1):end),Q_storage(plotParams,omega),'r-','linewidth',7)
        plot(f_hz_true(f_hz_true>=0),real(Q_hz_true(f_hz_true>=0)),'g-','linewidth',3)
        xlabel('Frequency [Hz]')
        ylabel('Storage Modulus [Pa]')
        legend('Data','Fit','Truth','location','best')
        grid on
        hold off

        nexttile
        scatter(f_hz(f_hz>=0),abs(imag(Q_hz(f_hz>=0))),'bo','linewidth',3)
        hold on
        plot(f_hz((end-numel(omega)+1):end),Q_loss(plotParams,omega),'r-','linewidth',7)
        plot(f_hz_true(f_hz_true>=0),abs(imag(Q_hz_true(f_hz_true>=0))),'g-','linewidth',3)
        xlabel('Frequency [Hz]')
        ylabel('Loss Modulus [Pa]')
        legend('Data','Fit','Truth','location','best')
        grid on
        hold off

        saveas(hfig,[path filesep 'NoiseStudy-' saveLabel '.fig'])
        print(hfig,[path filesep 'NoiseStudy-' saveLabel '.png'],'-dpng','-r300');
        
        % Save the results to our output structure
        temp = struct();
        
        temp.dataInClean = dataInTrue;
        temp.Q_hz_clean = Q_hz_true;
        temp.f_hz_clean = f_hz_true;
        
        temp.F = F;
        temp.h = h;
        temp.t = t;
        temp.z = z;
        temp.d = d;
                
        temp.dataIn = dataIn{j};
        temp.bestParams = bestParamsMap{j};
        temp.Q_hz = relaxanceMap{j};
        temp.f_hz = frequencyMap{j};
        temp.alphaInit = alphaMap{j};
        temp.omega = omega;
        
        switch j
            case 1
                outStruct(i_dir).WhiteGaussian = temp;
            case 2
                outStruct(i_dir).PinkNoise = temp;
            case 3
                outStruct(i_dir).WhiteNoise = temp;
            case 4
                outStruct(i_dir).BrownNoise = temp;
            case 5
                outStruct(i_dir).BlueNoise = temp;
            case 6
                outStruct(i_dir).PurpleNoise = temp;
            otherwise
                continue;
        end
        
    end

end

% Save output structure
save([simFileFolder filesep 'NoiseStudy-SNR' num2str(SNR) '-' date '.mat'],'-v7.3');

winopen(simFileFolder);
