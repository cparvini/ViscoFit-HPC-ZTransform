function sse = SSE_Maxwell_Map_zTransform(data,params,smoothOpt,windowsize,elasticSetting,fluidSetting)
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
    
    % Get Z-Transform Data from Curve
    [Q_hz,F_t,~,h_t,~,f_hz,alphaInit] = zTransformCurve(data,smoothOpt,windowsize);
    
    % Define Equations
    Q_hz = Q_hz(f_hz>=0);
    f_hz = f_hz(f_hz>=0);
    omega = (2.*pi.*f_hz);

    % Get functions
    n_terms = (numel(params)-2)/2;
    [Q_storage,Q_loss] = getZModels('maxwell',n_terms,alphaInit);
    
    % Calculate harmonics
    Qloss_hz = abs(imag(Q_hz));
    Qstorage_hz = real(Q_hz);
    
%     % Check the output!
%     try 
%         figure(hfig);
%         clf;
%         tiledlayout(1,3)
%     catch
%         hfig = figure;
%         tiledlayout(1,3)
%     end
% 
%     nexttile
%     plot(data{4},data{3},'b-','linewidth',3)
%     hold on
%     plot(h_t,F_t,'r-','linewidth',3)
%     xlabel('Indentation [m]')
%     ylabel('Force [N]')
%     if data{9}
%         legend('Data',['Thin-Sample, ' smoothOpt])
%     else
%         legend('Data',smoothOpt)
%     end
%     hold off
% 
%     nexttile
%     plot(f_hz,real(Q_hz),'b-','linewidth',3)
%     hold on
%     plot(f_hz,Q_storage(params,omega),'g-','linewidth',3)
%     xlabel('Frequency [Hz]')
%     ylabel('Storage Modulus [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     nexttile
%     plot(f_hz,abs(imag(Q_hz)),'b-','linewidth',3)
%     hold on
%     plot(f_hz,Q_loss(params,omega),'g-','linewidth',3)
%     xlabel('Frequency [Hz]')
%     ylabel('Loss Modulus [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     disp('pause');
    
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