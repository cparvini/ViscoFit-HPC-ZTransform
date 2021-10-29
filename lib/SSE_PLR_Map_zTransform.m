function sse = SSE_PLR_Map_zTransform(data,params,varargin)
    %SSE_PLR_Map_zTransform Calculate the SSE for the PLR model using the Z
    %Transform Method
    %   Calculate the Sum of Squared Errors for the Power Law
    %   Rheology Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_PLR function above. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_PLR function - a
    %   Sum of Squared Errors for that particular pixel.

    % Solve for z-transform quantities.
    % First, we need to find the correct alpha value and store
    % it for later. Create the initial guess first, using the
    % values from the force array.
    time = data{1};
    tipSize = data{5};
    nu = data{6};
    tipGeom = data{7};
    thinSample = data{9};
    h_finite = data{10};

    alphaInit = NaN(2,1);
    alphaInit(1) = (1/numel(data{1}))*(real(log(data{3}(1))) ...
        / real(log(data{3}(end))));
    alphaInit(2) = (1/numel(data{1}))*(real(log(data{4}(1))) ...
        / real(log(data{4}(end))));
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
        F_t = data{3}; % forces
        h_t = data{4}.^(betaParam); % indentations
        
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
                F_hz = modelfun(alphaInit,smooth(data{1},F_t,windowsize)');
                h_hz = modelfun(alphaInit,smooth(data{1},h_t,windowsize)');
                                
            case 'g-hz'
                % Gaussian smoothing in frequency
                F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowsize);
                h_hz = smoothdata(modelfun(alphaInit,h_t),'gaussian',windowsize);
                                
            case 'ma-hz'
                % Moving Average smoothing in frequency
                F_hz = smooth(data{1},modelfun(alphaInit,F_t),windowsize)';
                h_hz = smooth(data{1},modelfun(alphaInit,h_t),windowsize)';                
                
        end
        
        dt = mode(data{2});
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
%         scatter(data{1},data{3}./BEC,'bo','linewidth',3)
%         hold on
%         plot(data{1},F_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(data{1},data{4}.^(betaParam),'bo','linewidth',3)
%         hold on
%         plot(data{1},h_t,'r-','linewidth',3)
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
%         scatter(data{1},modelfun(alphaInit,data{3}./BEC),'bo','linewidth',3)
%         hold on
%         plot(data{1},F_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(data{1},modelfun(alphaInit,data{4}.^(betaParam)),'bo','linewidth',3)
%         hold on
%         plot(data{1},h_hz,'r-','linewidth',3)
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
%         scatter(f_hz,(modelfun(alphaInit,data{3}./BEC)./modelfun(alphaInit,data{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%         hold on
%         plot(f_hz,Q_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         xlim([min(f_hz) max(f_hz)])
%         grid on
%         hold off
%         legend('Original','Smooth','location','best')
        
    else
        
        h_t = data{4}.^(betaParam); % indentations
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
        
        F_t = data{3}./BEC; % forces
                        
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
                F_hz = modelfun(alphaInit,smooth(data{1},F_t,windowsize)');
                h_hz = modelfun(alphaInit,smooth(data{1},h_t,windowsize)');
                                
            case 'g-hz'
                % Gaussian smoothing in frequency
                F_hz = smoothdata(modelfun(alphaInit,F_t),'gaussian',windowsize);
                h_hz = smoothdata(modelfun(alphaInit,h_t),'gaussian',windowsize);
                                
            case 'ma-hz'
                % Moving Average smoothing in frequency
                F_hz = smooth(data{1},modelfun(alphaInit,F_t),windowsize)';
                h_hz = smooth(data{1},modelfun(alphaInit,h_t),windowsize)';                
                
        end
        
        dt = mode(data{2});
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
%         scatter(data{1},data{3}./BEC,'bo','linewidth',3)
%         hold on
%         plot(data{1},F_t,'r-','linewidth',3)
%         xlabel('Time [s]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(data{1},data{4}.^(betaParam),'bo','linewidth',3)
%         hold on
%         plot(data{1},h_t,'r-','linewidth',3)
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
%         scatter(data{1},modelfun(alphaInit,data{3}./BEC),'bo','linewidth',3)
%         hold on
%         plot(data{1},F_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Force [N]')
%         legend('Original','Smooth','location','best')
%         grid on
%         hold off
%         
%         nexttile
%         scatter(data{1},modelfun(alphaInit,data{4}.^(betaParam)),'bo','linewidth',3)
%         hold on
%         plot(data{1},h_hz,'r-','linewidth',3)
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
%         scatter(f_hz,(modelfun(alphaInit,data{3}./BEC)./modelfun(alphaInit,data{4}.^(betaParam))) .* (1./c),'bo','linewidth',3)
%         hold on
%         plot(f_hz,Q_hz,'r-','linewidth',3)
%         xlabel('Frequency [Hz]')
%         ylabel('Relaxance [Pa]')
%         xlim([min(f_hz) max(f_hz)])
%         grid on
%         hold off
%         legend('Original','Smooth','location','best')

    end
    
    % Calculate harmonics
    Qloss_hz = abs(imag(Q_hz));
    Qstorage_hz = real(Q_hz);
    
    % There is no explicit definition of "relaxance" for a PLR model in the
    % frequency domain in the publications we have found. Instead, we
    % calculate the time-domain modulus and transform that to the Z-domain
    % such that it can be compared to the storage and loss moduli from the
    % Z-transform of our data streams. It is crude, but functional.
    
    % Make our Dirac Delta array for the convolution portion
    diracArray = zeros(size(time));
    diracArray(time-dt<2*eps) = 1;

    % Set t' equal to dt, per Efremov (2017)
    t_prime = dt;

    % Calculate Power Law Rheology Modulus
    if length(params) > 1
        E_t = params(1).*( (1 + time./t_prime) .^ (-params(2)) );
    else
        E_t = params(1).*diracArray./dt;
    end
    
    E_hz = modelfun(alphaInit,E_t);
    
    

    % calculate global residual
    sse = sum((vertcat(Qstorage_hz,Qloss_hz)-vertcat(real(E_hz),abs(imag(E_hz)))).^2,'all');

    % Power Law Rheology Roster:
    % [E_0 alpha]
    ub = [1e12;1];
    lb = [1e-2;0];

    if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
        sse = Inf;
    end
    
end % End PLR SSE Z-Transform Map Function