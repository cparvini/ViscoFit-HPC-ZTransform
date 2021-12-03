function [Q_hz,F_t,F_hz,h_t,h_hz,f_hz,alphaInit] = zTransformCurve(dataIn,smoothOpt,windowsize,thinPixel)
%ZTRANSFORMCURVE Extract Viscoelastic Information from Input Data Using the
%Z-Transform Method on a Single Experiment
%   This function takes in the experiment data (dataIn), the smoothing
%   method to use for the observables, and the windowsize to use for the
%   smoothing.
%
%   dataIn: (cell array) AFM experimental data to use for the Z-Transform
%           method. The cell entries for "dataIn" should be (in order):
%           time, dt, force, indentation, tip size, nu (Poisson's ratio) of
%           the sample, tip geometry, minimum timescale for fitting, the
%           thin sample setting (T/F), and the sample thickness.
%
%   smoothOpt: (string) The smoothing operation to perform. Two options are
%           available---moving average, and gaussian. They can each be
%           performed in the time domain (pre-Z Transformation) or in the
%           modified Fourier domain (post-Z Transformation). The settings
%           are as follows: "ma-time", "ma-hz", "g-time", "g-hz", and of
%           course "none".
%
%   windowsize: (double) Size of the moving window to use for the
%           smoothing method. This is the fraction of the overall dataset
%           length to use for the moving average/gaussian smoothing. Must
%           be less than 50% (0.5).
%
%   thinPixel: (true/false) Whether to use thin sample correction on this
%           pixel. This is necessary for situations where we are not
%           ignoring the substrate, and we ARE using the thin-sample
%           correction for the cell. This prevents correcting the force
%           curve incorrectly.
%

if windowsize > 0.5
    error('Your specified "windowsize" is too large a fraction of the dataset. Please ensure it is less than 50% (i.e. 0.5).');
end

% Make sure we're handling rows
for i = 1:6
    if ~isrow(dataIn{i})
        dataIn{i} = dataIn{i}';
    end
end

% Solve for z-transform quantities.
% First, we need to find the correct alpha value and store
% it for later. Create the initial guess first, using the
% values from the force array.
tipSize = dataIn{5};
nu = dataIn{6};
tipGeom = dataIn{7};
thinSample = thinPixel;
h_finite = dataIn{10};
windowBins = floor((numel(dataIn{3})*windowsize)/2); % divide by 2 b/c we will only fit positive frequencies later
if windowBins < 1
    windowBins = 1;
end
    
alphaInit = NaN(2,1);
alphaInit(1) = (1/numel(dataIn{1}))*((log(dataIn{3}(1))) ...
    / (log(dataIn{3}(end))));
alphaInit(2) = (1/numel(dataIn{1}))*(real(log(dataIn{4}(1))) ...
    / real(log(dataIn{4}(end))));
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
    case "punch"
        % Lopez-Guerra & Solares (Beilstein J. of Nanotech., 2017)
        c = (4*tipSize)./(1-nu);
        betaParam = 1;
    case "pyramidal"
        % Weber, Iturri, Benitez, and Toca-Herrera (Microscopy Research and Technique, 2019)
        % Modified Sneddon
        c = (tan(tipSize.*pi./180))./(sqrt(2).*(1-nu.^2));
        betaParam = 2;
end

% Optimize alpha
modelfun = @(alpha,x) (1/numel(x)).*fftshift(fft(x.*(exp(-alpha.*(1:numel(x))))));
F_hz = NaN;
h_hz = NaN;
if ~thinSample
    F_t = dataIn{3}; % forces
    h_t = dataIn{4}; % indentations

    switch smoothOpt
        case 'none'
            % Do not smooth
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'g-time'
            % Gaussian smoothing in time
            F_t = smoothdata(F_t,'gaussian',windowBins)';
            h_t = smoothdata(h_t,'gaussian',windowBins)';
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'ma-time'
            % Moving Average smoothing in time
            F_t = smooth(F_t,windowBins)';
            h_t = smooth(h_t,windowBins)';
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'g-hz'
            % Gaussian smoothing in frequency
            F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowBins);
            h_hz = smoothdata(modelfun(alphaInit,h_t.^betaParam),'gaussian',windowBins);

        case 'ma-hz'
            % Moving Average smoothing in frequency
            F_hz = smooth(modelfun(alphaInit,F_t),windowBins)';
            h_hz = smooth(modelfun(alphaInit,h_t.^betaParam),windowBins)';                

    end

    dt = mode(dataIn{2});
    fs = dt.^-1;
    N = length(F_hz);
    df = fs/N;
    f_hz = -fs/2:df:fs/2-df + (df/2)*mod(N,2);
    f_hz(abs(f_hz) < 1e-6) = 0;

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
%         scatter(dataIn{1},dataIn{3}./BEC,'bo','linewidth',3)
%         hold on
%         plot(dataIn{1},F_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{1},dataIn{4},'bo','linewidth',3)
%         hold on
%         plot(dataIn{1},h_t,'r-','linewidth',3)
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
%         scatter(dataIn{1},modelfun(alphaInit,dataIn{3}./BEC),'bo','linewidth',3)
%         hold on
%         plot(dataIn{1},F_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(dataIn{1},modelfun(alphaInit,dataIn{4}.^(betaParam)),'bo','linewidth',3)
%         hold on
%         plot(dataIn{1},h_hz,'r-','linewidth',3)
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
%         scatter(f_hz,(modelfun(alphaInit,dataIn{3}./BEC)./modelfun(alphaInit,dataIn{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%         hold on
%         plot(f_hz,Q_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         xlim([min(f_hz) max(f_hz)])
%         grid on
%         hold off
%         legend('Original','Smooth','location','best')

else

    h_t = dataIn{4}; % indentations
    BEC = NaN;

    switch tipGeom
        case "spherical"
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            a = sqrt(tipSize.*h_t);
            BEC = (1/(h_finite^0)) + ((1.133.*(a))./(h_finite)) + ...
                ((1.497.*(a.^2))./(h_finite.^2)) + ...
                ((1.469.*(a.^3))./(h_finite.^3)) + ...
                ((0.755.*(a.^4))./(h_finite.^4));

        case "conical"
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            a = tan(tipSize.*pi./180);
            BEC = (1/(h_finite^0)) + ((0.721.*(h_t).*(a))./(h_finite)) + ...
                ((0.650.*(h_t.^2).*(a.^2))./(h_finite^2)) + ...
                ((0.491.*(h_t.^3).*(a.^3))./(h_finite^3)) + ...
                ((0.225.*(h_t.^4).*(a.^4))./(h_finite^4));

        case "punch"
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            % Note: this correction does not depend on the indentation
            % depth, according to the Garcia & Garcia paper (Results).
            a = tipSize;                                                % Radius of the punch
            BEC = (1/(h_finite^0)) + ((1.133*(a))/(h_finite)) + ...
                ((1.283*(a^2))/(h_finite^2)) + ...
                ((0.598*(a^3))/(h_finite^3)) - ...
                ((0.291*(a^4))/(h_finite^4));

        case "pyramidal"
            % ASSUMPTION: BEC is the same for Conical and Pyramidal Probes
            % If one can find a pyramidal BEC, replace this
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            a = tan(tipSize.*pi./180);
            BEC = (1/(h_finite^0)) + ((0.721.*(h_t).*(a))./(h_finite)) + ...
                ((0.650.*(h_t.^2).*(a.^2))./(h_finite^2)) + ...
                ((0.491.*(h_t.^3).*(a.^3))./(h_finite^3)) + ...
                ((0.225.*(h_t.^4).*(a.^4))./(h_finite^4));

    end

    F_t = dataIn{3}./BEC; % forces

    % Perform Smoothing
    
    switch smoothOpt
        case 'none'
            % Do not smooth
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'g-time'
            % Gaussian smoothing in time
            F_t = smoothdata(F_t,'gaussian',windowBins)';
            h_t = smoothdata(h_t,'gaussian',windowBins)';
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'ma-time'
            % Moving Average smoothing in time
            F_t = smooth(F_t,windowBins)';
            h_t = smooth(h_t,windowBins)';
            F_hz = modelfun(alphaInit,F_t);
            h_hz = modelfun(alphaInit,h_t.^betaParam);

        case 'g-hz'
            % Gaussian smoothing in frequency
            F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowBins);
            h_hz = smoothdata(modelfun(alphaInit,h_t.^betaParam),'gaussian',windowBins);

        case 'ma-hz'
            % Moving Average smoothing in frequency
            F_hz = smooth(modelfun(alphaInit,F_t),windowBins)';
            h_hz = smooth(modelfun(alphaInit,h_t.^betaParam),windowBins)';                

    end

    dt = mode(dataIn{2});
    fs = dt.^-1;
    N = length(F_hz);
    df = fs/N;
    f_hz = -fs/2:df:fs/2-df + (df/2)*mod(N,2);
    f_hz(abs(f_hz) < 1e-6) = 0;

    Q_hz = (F_hz./h_hz) .* (1./c);

%     % Observe Smoothing Results!
%     try
%         figure(hfig)
%         clf;
%     catch
%         hfig = figure;
%         tiledlayout(2,1);
%     end
%     nexttile
%     scatter(dataIn{1},dataIn{3}./BEC,'bo','linewidth',3)
%     hold on
%     plot(dataIn{1},F_t,'r-','linewidth',3)
%     xlabel('Time [s]')
%     ylabel('Force [N]')
%     legend('Original','Smooth','location','best')
%     grid on
%     hold off
% 
%     nexttile
%     scatter(dataIn{1},dataIn{4},'bo','linewidth',3)
%     hold on
%     plot(dataIn{1},h_t,'r-','linewidth',3)
%     xlabel('Time [s]')
%     ylabel('Indentation [m]')
%     legend('Original','Smooth','location','best')
%     grid on
%     hold off
% 
%     try
%         figure(hzfig)
%         clf;
%     catch
%         hzfig = figure;
%         tiledlayout(2,1);
%     end
%     nexttile
%     scatter(f_hz,modelfun(alphaInit,dataIn{3}./BEC),'bo','linewidth',3)
%     hold on
%     plot(f_hz,F_hz,'r-','linewidth',3)
%     xlabel('Frequency [Hz]')
%     ylabel('Force [N]')
%     legend('Original','Smooth','location','best')
%     grid on
%     hold off
% 
%     nexttile
%     scatter(f_hz,modelfun(alphaInit,dataIn{4}.^(betaParam)),'bo','linewidth',3)
%     hold on
%     plot(f_hz,h_hz,'r-','linewidth',3)
%     xlabel('Frequency [Hz]')
%     ylabel('Indentation [m]')
%     legend('Original','Smooth','location','best')
%     grid on
%     hold off
% 
%     try
%         figure(Qhzfig)
%         clf;
%     catch
%         Qhzfig = figure;
%     end
%     scatter(f_hz,(modelfun(alphaInit,dataIn{3}./BEC)./modelfun(alphaInit,dataIn{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%     hold on
%     plot(f_hz,Q_hz,'r-','linewidth',3)
%     xlabel('Frequency [Hz]')
%     ylabel('Relaxance [Pa]')
%     xlim([min(f_hz) max(f_hz)])
%     grid on
%     hold off
%     legend('Original','Smooth','location','best')

end

end

