function [] = createTestMap(startDir,saveDir,saveLabel,varargin)
%CREATETESTMAP Create A Simulated QI Map Using Simulated SFS Results
%   Taking in a directory (startDir), this function will create a
%   bullseye-shaped map using all of the simulation results available
%   inside of startDir. Optionally, if the user wishes, they can pass a
%   noise color and magnitude as additional arguments. Importantly, unless
%   otherwise specified as the first optional argument, the map size will
%   default to 128x128 pixels. All files will be saved with the "saveLabel"
%   string added to their filename.

mapSize = [128 128];    % Output Map Size
noiseType = 'none';     % Noise Color
noiseMag = 50;          % SNR
if nargin > 3
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        mapSize = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        noiseType = varargin{i};                        
                    end
                case 3
                    if ~isempty(varargin{i})
                        noiseMag = varargin{i};                        
                    end
                otherwise
                    fprintf('Passed additional parameters to createTestMap() which were not used.');
            end
        end
    end
end

% Analysis Settings
includeRetract = 0;             % Include data from the retract curve
filterType = 'none';            % Choose the filter used to smooth data
N = 2;                          % Order of Butterworth filter, if used
cutoff_Hz = 5000;               % Cutoff frequency
findRep = 'legendre';           % Search direction for the repulsive region
maxHeight = 10e-6;              % When making the fake height map, what is the max height?

fprintf('Creating Simulated Test Map...\n');

Files=dir([startDir filesep '*.*']);

% Remove Directories
Files=Files(~ismember({Files.name},{'.','..'}));
toRemove = find([Files.isdir] == 1);
Files(toRemove) = [];

% Remove Filetypes To Ignore
toRemove = find(~endsWith({Files.name}, {'.mat'}));
Files(toRemove) = [];

toRemove = find(contains({Files.name}, {'settingsStruct','Settings'}));
Files(toRemove) = [];

toRemove = find(contains({Files.name}, {'FitResults','MapResults','PlotResults','mapStruct'}));
Files(toRemove) = [];

if isempty(Files)
    error('ERROR: No valid AFM files found! Please ensure you have selected the correct directory, and that LoadAFMData contains a case describing how to load that filetype.');
end

fileLabels = {};
for i = 1:length(Files)
    FilesTempPrep = Files(i).name;
    FilesTempPrep = strsplit(FilesTempPrep, {'_' '.' ' '},'CollapseDelimiters',true);
    fileLabels = FilesTempPrep{1};
end

% Loop through the files and proceed with loading, parsing, and calculation
% for each one.
dataStruct = struct;
for k = 1:length(Files)
    
    FileNames = Files(k).name;
    FileInfo = strsplit(FileNames, {'_' '.' ' '},'CollapseDelimiters',true);
    
    strSwitch = FileInfo{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end
    
    % This case is for loading the "TestCondition" simulation files
    % which were generated for testing the code bases. The files
    % have a regular, relatively simple format and have two
    % separate .mat files for each: one contains the data from the
    % simulation (loaded here into RawData), and the other is a
    % settings file (used to perform the simulation, stored here in
    % settingsData).
    RawData = load([Files(k).folder filesep Files(k).name],'z','d','t','F');

    settingsCheck=dir([startDir '/*.*']);
    
    % Remove Directories
    settingsCheck=settingsCheck(~ismember({settingsCheck.name},{'.','..'}));
    toRemove = find([settingsCheck.isdir] == 1);
    settingsCheck(toRemove) = [];

    % Remove Filetypes To Ignore
    toRemove = find(~endsWith({settingsCheck.name}, {'.mat'}));
    settingsCheck(toRemove) = [];

    toRemove = find(~contains({settingsCheck.name}, {'Settings'}));
    settingsCheck(toRemove) = [];

    if size(settingsCheck,1) > 1
        setInd = k;
    else
        setInd = 1;
    end

    settingsData = load([settingsCheck(setInd).folder filesep settingsCheck(setInd).name]);
    dataStruct(k).tipGeom = settingsData.settingsStruct.tipGeom;
    dataStruct(k).z = -(RawData.z-RawData.z(1));
    dataStruct(k).t = RawData.t;
    dataStruct(k).n_sam = numel(dataStruct(k).t);

    k_cantilever(k) = settingsData.settingsStruct.k_m1;
    dataStruct(k).k_cantilever = k_cantilever(k);
    dataStruct(k).r_tip = settingsData.settingsStruct.r_tip;
    dataStruct(k).nu_sample = settingsData.settingsStruct.nu_sample;
    dataStruct(k).d = RawData.F./k_cantilever(k);
    defl_InVOLS(k) = 1;
    dataStruct(k).dt = settingsData.settingsStruct.dt;
    
    dataStruct(k).trueParams = [settingsData.settingsStruct.E0;settingsData.settingsStruct.phi_f];
    trueParams = [];
    for j = 1:size(settingsData.settingsStruct.linearParamData,1)
        trueParams = vertcat(trueParams,...
            [settingsData.settingsStruct.linearParamData(j,2);...
            settingsData.settingsStruct.linearParamData(j,4)]);
    end
    
    dataStruct(k).thinSample = false;

%     % Quick Look
%     plot(dataStruct(k).t,RawData.d)
    
end

% Analyze the results we have extracted from the simulation files
for k = 1:size(dataStruct,2)
        
    % Find the approach portion of the data    
    [~, z_max_ind] = max(dataStruct(k).z);
    
    % Make sure we are handling row vectors
    if ~isrow(dataStruct(k).z) dataStruct(k).z = dataStruct(k).z'; end
    if ~isrow(dataStruct(k).d) dataStruct(k).d = dataStruct(k).d'; end
    if ~isrow(dataStruct(k).t) dataStruct(k).t = dataStruct(k).t'; end
    
    % Trim the data. If we include the retract portion of the dataset, then
    % we will go until the applied force detected is negative. If it is not
    % included, we end at the point of maximum z-sensor (maximum
    % deformation condition).
    if ~includeRetract
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    else
        F_temp = dataStruct(k).d(z_max_ind:end) .* k_cantilever(k);
        if ~isempty(F_temp) && length(F_temp) > 1
            non_contact_ind = find(F_temp < 0,1);
            if isempty(non_contact_ind) non_contact_ind = length(F_temp)-1; end
        else
            non_contact_ind = 0;
        end
        z_max_ind = non_contact_ind+z_max_ind;
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    end
    
    % Filter and Shift z_sensor data
    switch filterType
        case 'butter'
            % Create the butterworth
            [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(k).dt)),'low'); % This makes a lowpass filter

            d_smooth = (filter(b,a,dataStruct(k).d_approach));             % Next, apply the filter
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end
            
            % Calculate Deflection Offset
            [~, d_min_ind] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
        case {'FIR','IIR'}
            
            Fs = 1/(dataStruct(k).dt); 
            Fstop = ( 1 / dataStruct(k).dt )/10;
            if Fstop >= 1/(2*dataStruct(k).dt)
                Fstop = 1/(2.05*dataStruct(k).dt);
            elseif Fstop < 1/(10*dataStruct(k).dt)
                Fstop = 1/(10*dataStruct(k).dt);
            end
            Fpass = Fstop*0.01;
            Rp = 0.01;
            Astop = 80;

            Nfilter = round((Fs/(Fstop-Fpass))*(Astop/22)); % "Fred Harris Rule of Thumb" for Filter Order
            firlo = fir1(Nfilter,0.5,'low');
            d_smooth = filter(firlo,1,dataStruct(k).d_approach);
            delay = floor(mean(grpdelay(firlo,Nfilter)));
            
            sf = d_smooth;
            if delay < length(sf)
                % Correct filter delay
                sf(1:delay) = [];
            else
                delay = 0;
            end

            d_smooth = sf;
            
            % Calculate Deflection Offset
            [~, dSmoothMin] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(dSmoothMin*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected((delay+1):end);

        case 'none'
            d_smooth = dataStruct(k).d_approach;
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end
            
            % Calculate Deflection Offset
            [~, d_min_ind] = min(d_smooth);
            indScale = 0.9;
            d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
            dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
            dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
    end
    
    % Correct the z-sensor data
    dataStruct(k).z_full_corrected = dataStruct(k).z - ...
        dataStruct(k).z_approach(dSmoothMin) + ...
        dataStruct(k).d_corrected(dSmoothMin);

    % Store our detected offsets
    dataStruct(k).d0 = dataStruct(k).d_corrected(dSmoothMin);
    dataStruct(k).z0 = dataStruct(k).z_corrected(dSmoothMin);
    
    % Calculate Force and Indentation
    dataStruct(k).F = dataStruct(k).d_corrected .* k_cantilever(k);
    dataStruct(k).F_smooth = d_smooth .* k_cantilever(k);
    dataStruct(k).d_smooth = d_smooth;
    dataStruct(k).z_smooth = z_smooth;
    dataStruct(k).h = (dataStruct(k).z(1:z_max_ind) - dataStruct(k).z0)...
        - (dataStruct(k).d(1:z_max_ind) - dataStruct(k).d0);
    
    % Get Repulsive Force Application Portion of the Data
    n_offset = length(dataStruct(k).d_corrected(dSmoothMin:z_max_ind));
    n_offset_smooth = length(dataStruct(k).d_corrected(dSmoothMin:(z_max_ind-delay)));

    if n_offset < 10
        fprintf('\nBad Dataset: %s\nStoring fake data here.\n',fileLabels{k});
        
        % Store fake arrays so the script doesn't fail for this bad pixel.
        dataStruct(k).t_r = dataStruct(k).t;
        dataStruct(k).z_r = NaN(size(dataStruct(k).z));
        dataStruct(k).d_r = NaN(size(dataStruct(k).d));
        dataStruct(k).t_r_smooth = dataStruct(k).t_r;
        dataStruct(k).z_r_smooth = NaN(size(dataStruct(k).z));
        dataStruct(k).d_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).tip_rep_pos = 1;
        dataStruct(k).tip_rep_pos_smooth = 1;
        tip_rep_pos_all(k) = 1;
        dSmoothMinAll(k) = 1;
        dataStruct(k).F_r = NaN(size(dataStruct(k).d));
        dataStruct(k).h_r = NaN(size(dataStruct(k).d));
        dataStruct(k).F_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).h_r_smooth = NaN(size(dataStruct(k).d));
        dataStruct(k).z_max_ind = length(dataStruct(k).z);
        dataStruct(k).z_max_ind_smooth = length(dataStruct(k).z);
        dataStruct(k).dSmoothMin = 1;
        dataStruct(k).F_r_log = NaN(size(dataStruct(k).d));
        dataStruct(k).t_r_log = NaN(size(dataStruct(k).d));
        
        continue;

    end
    
    dt = dataStruct(k).dt;    
    
    % Create a "clean" repulsive time array
    t_rep = linspace(0,(n_offset-1)*dt,n_offset);
    t_rep_smooth = linspace(0,(n_offset_smooth-1)*dt,n_offset_smooth);

    % Store the repulsive z-sensor and deflection for normal and smooth
    % cases.
    z_rep = dataStruct(k).z_corrected(dSmoothMin:end);
    d_rep = dataStruct(k).d_corrected(dSmoothMin:end);
    z_rep_smooth = dataStruct(k).z_smooth(dSmoothMin:end);
    d_rep_smooth = dataStruct(k).d_smooth(dSmoothMin:end);
    
    % Store the tip position
    tip_rep = d_rep;
    tip_rep_smooth = d_rep_smooth;
    
    % Search for the repulsive portion of the data starting from either the
    % first index and moving to the end ('forward') or starting at the end
    % and walking backward ('reverse') and looking for where the data is
    % appreciably above the approach data. Alternately, perform a Legendre
    % transformation on the data with a moving window (modified Wavelet
    % Transform) to remove noise and make the repulsive regime obvious.
    if strcmp(findRep,'legendre')
    
        % Use the Legendre-Wavelet approach
        data = d_rep.*k_cantilever(k);
        if ~includeRetract
            % Make the data symmetric
            data = horzcat(data, flip(data));
        end
        
        n_data = numel(data);
        n_legs = 10;
        n_window = 10;

        data_split = NaN((n_data-n_window),n_window);
        ai = 0;
        bi = n_window-1;
        for ii = 1:size(data_split,1)
            ai = ai+1;
            bi = bi+1;
            if bi < numel(data)
                data_split(ii,:) = data(ai:bi);
            else
                break;
            end
        end

        L_out = NaN(n_window,n_legs);
        x_temp = linspace(-1,1,size(data_split,2));
        L_func = legendre(n_legs-1,x_temp);
        L_smooth = NaN(size(L_out,1)-1,1);

        for ii = 1:size(data_split,1)
            for jj = 1:n_legs
                % Legendre Transform
                L_out(ii,jj) = dot(data_split(ii,:),L_func(jj,:));
            end
            if ii > 1
                L_smooth(ii-1) = dot(L_out(ii,:),L_out(ii-1,:));
            end
        end

        data = (data(1:numel(data)/2));
        L_smooth = (L_smooth(1:ceil(numel(L_smooth)/2)));

        try
            temp = 5*max(abs(L_smooth(1:ceil(numel(L_smooth)/20))));
            upperLim = mean(L_smooth(1:ceil(numel(L_smooth)/20)))+temp;
            idx = find(L_smooth>upperLim,1,'first')-1;
        catch
            idx = [];
        end
        
        % Prevent removing too much data
        tip_rep_pos = max([idx 2]);
        tip_rep_pos_smooth = max([idx-delay 2]);

%         % Visualize the contact point estimate!
%         try 
%             figure(tfig);
%             clf
%             tiledlayout(2,1);
%         catch
%             tfig = figure;
%             tiledlayout(2,1);
%         end
% 
%         nexttile
%         plot(1:numel(L_smooth),L_smooth,'r-','linewidth',3)
%         hold on
%         scatter(idx,L_smooth(idx),'bo','linewidth',3)
%         plot(1:numel(L_smooth),upperLim*ones(size(1:numel(L_smooth))),'-b','linewidth',3)
%         xlabel('Bin Number')
%         ylabel('Magnitude')
%         xlim([1 numel(L_smooth)])
%         title('Smooth Bin Data')
%         hold off
% 
%         nexttile
%         plot(1:numel(data),data,'r-','linewidth',3)
%         hold on
%         scatter(idx,data(idx),'bo','linewidth',3)
%         xlabel('Data Index')
%         ylabel('Force [N]')
%         xlim([1 numel(data)])
%         title('Original Force Data')
%         hold off
% 
%         disp('pause');
        
    else
    
        % Legacy Method
        if strcmp(findRep,'forward')
            tip_rep_pos = find(tip_rep>0,1);                                   % Find first position above 0
            if isempty(tip_rep_pos) || tip_rep_pos < 2
                tip_rep_pos = 2;
            end
            
            tip_rep_pos_smooth = find(tip_rep_smooth>0,1);                     % Find first position above 0
            if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth < 2
                tip_rep_pos_smooth = 2;
            end
        elseif strcmp(findRep,'reverse')
            tip_rep_pos = (length(tip_rep) - find(flip(tip_rep)<0,1));       % Find last position above 0
            if isempty(tip_rep_pos) || tip_rep_pos < 2
                tip_rep_pos = 2;
            end
            
            tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flip(tip_rep_smooth)<0,1));   % Find last position above 0
            if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth < 2
                tip_rep_pos_smooth = 2;
            end
        end
        
    end
    
    % Store the repulsive force application portion of the dataset.
    dataStruct(k).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos-1);
    dataStruct(k).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos-1);
    dataStruct(k).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos-1);
    dataStruct(k).t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth-1);
    dataStruct(k).z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth-1);
    dataStruct(k).d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth-1);
    dataStruct(k).tip_rep_pos = tip_rep_pos;
    dataStruct(k).tip_rep_pos_smooth = tip_rep_pos_smooth;
    
    % Store for benchmarking later. This will help us decide which region
    % of the curves we can actually "average" together, since we want to
    % use data where all files that are averaged have some contribution
    % (i.e. where they overlap).
    tip_rep_pos_all(k) = tip_rep_pos;
    dSmoothMinAll(k) = dSmoothMin;
    
    % Calculate Force and Indentation during Repulsive Portion
    dataStruct(k).F_r = dataStruct(k).d_r .* k_cantilever(k); % Calculate Force
    dataStruct(k).h_r = (dataStruct(k).z_r - dataStruct(k).d_r); % Calculate Indentation
    
    dataStruct(k).F_r_smooth = dataStruct(k).d_r_smooth .* k_cantilever(k); % Calculate Smooth Force
    dataStruct(k).h_r_smooth = dataStruct(k).z_r_smooth - dataStruct(k).d_r_smooth; % Calculate Smooth Indentation
    
    % Store our critical indices
    dataStruct(k).z_max_ind = z_max_ind;
    dataStruct(k).z_max_ind_smooth = z_max_ind-delay;
    dataStruct(k).dSmoothMin = dSmoothMin;

    % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
    tr = dataStruct(k).dt;
    st = dataStruct(k).t(z_max_ind);
    
    F_log = [];
    t_log = [];
    
    F_log = log_scale(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
    t_log = log_scale(dataStruct(k).t_r,dataStruct(k).t_r,tr,st);
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(k).F_r_log = F_log;
    dataStruct(k).t_r_log = t_log;

end

% Begin creating the map datasets that we will need for our cluster script
n_pixels = prod(mapSize);
forces = cell(1,n_pixels);
times = cell(1,n_pixels);
indentations = cell(1,n_pixels);
tipSize = cell(1,n_pixels);
nu = cell(1,n_pixels);
mapSizeAll = cell(1,n_pixels);
pixelHeight = cell(1,n_pixels);
trueParams = cell(1,n_pixels);

% Create a roster of pixels
xdata = 1:mapSize(1);
ydata = flip(1:mapSize(2));
[X, Y] = meshgrid(xdata,ydata);
X = X - mapSize(1)/2;
Y = Y - mapSize(2)/2;
distMat = sqrt(X.^2 + Y.^2);
heightMat = abs(distMat-max(distMat,[],'all')) ./ max(distMat,[],'all');
sections = linspace(0,1,1+size(dataStruct,2)); % Create our binning criteria

binRoster = ones(mapSize);
for k = 1:size(dataStruct,2)
    ai = sections(k)*max(distMat,[],'all');
    bi = sections(k+1)*max(distMat,[],'all');
    if k < size(dataStruct,2)
        binRoster(distMat >= ai & distMat < bi) = numel(sections)-k;
    else
        binRoster(distMat >= ai & distMat <= bi) = numel(sections)-k;
    end
end

xc = 1;
yc = 0;
for k = 1:n_pixels
        
    % Get the current pixel position
    if xc > mapSize(1)
        xc = 1;
        yc = yc + 1;
    end
    idx_pixel = sub2ind(flip(mapSize),mapSize(2)-yc,xc);
    rowID = binRoster(idx_pixel);
    
    % Put the raw data in the right place
    tipSize{k} = dataStruct(rowID).r_tip;
    nu{k} = dataStruct(rowID).nu_sample;
    mapSizeAll{k} = mapSize;
    pixelHeight{k} = heightMat(idx_pixel).*maxHeight;
    trueParams{k} = dataStruct(rowID).trueParams;
    
    % Settings
    SNR = noiseMag;
    signalPower = sum(dataStruct(rowID).F_r.^2)/numel(dataStruct(rowID).F_r);
    
    switch noiseType
        case 'awgn'
            % white gaussian noise
            forces{k} = awgn(dataStruct(rowID).F_r,SNR,'measured');
            
        case 'pink'
            temp = dsp.ColoredNoise('pink',numel(dataStruct(rowID).F_r),1);
            noise = temp();
            if ~isrow(noise)
                noise = noise';
            end
            noisePower = sum(noise.^2)/numel(noise);
            scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
            forces{k} = dataStruct(rowID).F_r + noise*scaleFactor;
            
        case 'white'
            temp = dsp.ColoredNoise('white',numel(dataStruct(rowID).F_r),1);
            noise = temp();
            if ~isrow(noise)
                noise = noise';
            end
            noisePower = sum(noise.^2)/numel(noise);
            scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
            forces{k} = dataStruct(rowID).F_r + noise*scaleFactor;
            
        case 'brown'
            temp = dsp.ColoredNoise('brown',numel(dataStruct(rowID).F_r),1);
            noise = temp();
            if ~isrow(noise)
                noise = noise';
            end
            noisePower = sum(noise.^2)/numel(noise);
            scaleFactor = sqrt(signalPower./(noisePower*(10^(SNR/10))));
            forces{k} = dataStruct(rowID).F_r + noise*scaleFactor;
            
        case 'blue'
            temp = dsp.ColoredNoise('blue',numel(dataStruct(rowID).F_r),1);
            noise = temp();
            if ~isrow(noise)
                noise = noise';
            end
            noisePower = sum(noise.^2)/numel(noise);
            scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
            forces{k} = dataStruct(rowID).F_r + noise*scaleFactor;
            
        case 'purple'
            temp = dsp.ColoredNoise('purple',numel(dataStruct(rowID).F_r),1);
            noise = temp();
            if ~isrow(noise)
                noise = noise';
            end
            noisePower = sum(noise.^2)/numel(noise);
            scaleFactor = sqrt(signalPower/(noisePower*(10^(SNR/10))));
            forces{k} = dataStruct(rowID).F_r + noise*scaleFactor;
            
        otherwise
            forces{k} = dataStruct(rowID).F_r;
    end
    
    times{k} = dataStruct(rowID).t_r;
    indentations{k} = dataStruct(rowID).h_r;
    
    xc = xc + 1;
    
end

% Start Timer
tic

% Make a structure for our settings
classSettings = struct;

classSettings.nu = nu;                                      % Sample Poisson Ratio for all curves
classSettings.tipGeom = dataStruct(1).tipGeom;              % Tip geometry for these experiments
classSettings.thinSample = dataStruct(1).thinSample;        % Use the thin-sample correction feature
classSettings.correctTilt = false;                          % Correct any tilt in the substrate by fitting a plane to the image edges. Do NOT use for monolayers.
classSettings.hideSubstrate = false;                        % Remove pixels designated as "substrate" from the visco fitting
classSettings.zeroSubstrate = true;                         % In addition to correcting for tilt, also set the new, flat surface to have a minimum value starting at zero.
classSettings.optimizeFlattening = true;                    % Attempt to search for the best-order polynomial surface for correction

classSettings.pixelHeight = pixelHeight;                    % The height for each pixel in the map. Used for thresholding/thin sample correction

% Create the class object
viscoZ = ViscoFitZ(forces,times,indentations,tipSize,classSettings);

% Make a structure for our settings
extractionSettings = struct;
extractionSettings.N_workers = 2;                           % Pass the number of workers to the fitting
extractionSettings.hideSubstrate = viscoZ.hideSubstrate;    % Remove pixels designated as "substrate" from the visco fitting
extractionSettings.windowsize = 0.05;                       % Window size for smoothing methods
extractionSettings.smoothOpt = 'none';                      % Which smoothing setting to use on the harmonics.
                                                            % Options: none, g-time, ma-time, g-hz, ma-hz
                                                            % *-time smooths before z-transform. *-hz will
                                                            % smooth after z-transforming F and h.

zTransformAnalysis = processMapZ_func(viscoZ,extractionSettings);

% In case the class was unsuccessfully saved in
% processMapZ_func
if isempty(zTransformAnalysis.ViscoClass)
    fprintf('\n\nprocessMapZ_func resulted in an empty ViscoClass!\n\n');
    zTransformAnalysis.ViscoClass = viscoZ;
end

zTransformAnalysis.trueParamsMap = trueParams;
zTransformAnalysis.trueBinsMap = binRoster;
zTransformAnalysis.mapSize = mapSize;
zTransformAnalysis.windowsize = extractionSettings.windowsize;
zTransformAnalysis.smoothOpt = extractionSettings.smoothOpt;
zTransformAnalysis.correctTilt = classSettings.correctTilt;
zTransformAnalysis.hideSubstrate = extractionSettings.hideSubstrate;
zTransformAnalysis.thinSample = classSettings.thinSample;
save([saveDir filesep saveLabel '_MapResults_zTransform.mat'],'zTransformAnalysis','-v7.3')

analysisTime = toc;
fprintf('\nZ-Transform Analysis Time: %4.2f Minutes\n',analysisTime/60);
fprintf('\nTest map data generation complete.\n\n');

end