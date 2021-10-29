function sse = SSE_Maxwell_Map_zTransform(data,params,elasticSetting,fluidSetting)
    %SSE_Maxwell Calculate the SSE for the Maxwell model using the
    %Z-Transform approach
    %   Calculate the Sum of Squared Errors for the Generalized
    %   Maxwell Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_Maxwell function. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_Maxwell function - a
    %   Sum of Squared Errors for that particular pixel.
    
    % Solve for z-transform quantities.
    % First, we need to find the correct alpha value and store
    % it for later. Create the initial guess first, using the
    % values from the force array.
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
    
    % Define Equations
    omega = (2.*pi.*f_hz);
%     z_model = alphaInit + 1i.*omega;
    if n_terms > 0
%         Q_model = @(c,z) (sum(c(1:2:end))-sum( (repmat(ones(size(z)),n_terms,1).*c(3:2:end).*(c(4:2:end)./dt))./((c(4:2:end)./dt).*(repmat(1-z.^-1,n_terms,1))+1) .* repmat(1-z.^-1,n_terms,1) , 1));
        Q_storage = @(c,omega) sum(c(1:2:end)) - sum( (c(3:2:end).*repmat(ones(size(omega)),n_terms,1).*(((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)))./(((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).*((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1))), 1);
        Q_loss = @(c,omega) sum( repmat(omega,n_terms,1).*(c(3:2:end).*(c(4:2:end)./dt).*repmat(ones(size(omega)),n_terms,1).*(((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)))./((((1+(c(4:2:end)./dt).*alphaInit).^2).*repmat(ones(size(omega)),n_terms,1)).*((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1))), 1);
    else
%         Q_model = @(c,z) c(1).*ones(size(omega));
        Q_storage = @(c,omega) c(1).*ones(size(omega));
        Q_loss = @(c,omega) zeros(size(omega));
    end
    
    % Calculate harmonics
    Qloss_hz = abs(imag(Q_hz));
    Qstorage_hz = real(Q_hz);
    
    % calculate global residual
    sse = sum((vertcat(Qstorage_hz,Qloss_hz)-vertcat(Q_storage(params,omega),Q_loss(params,omega))).^2,'all');

    [tauInds,modulusInds] = getParamIndices(params);
    ub = zeros(size(params))+eps;
    lb = zeros(size(params));

    ub(modulusInds) = 1e12;
    lb(modulusInds) = 1e-2;

    tauCenters = data{8}.*(10.^( (1:length(params(3:2:end)))-1 ));
    ub(tauInds) = tauCenters*10;
    lb(tauInds) = tauCenters/10;

    if length(params) > 1
        if fluidSetting
            ub(2) = max(tauCenters)*1e2;
            lb(2) = min(data{2});
        end
    end

    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
    
end % End Maxwell SSE Z-Transform Map Function