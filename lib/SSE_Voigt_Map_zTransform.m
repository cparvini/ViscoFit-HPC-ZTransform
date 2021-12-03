function sse = SSE_Voigt_Map_zTransform(data,params,ub,lb,smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting)
    %SSE_Voigt Calculate the SSE for the Voigt model using the Z-Transform
    %approach
    %   Calculate the Sum of Squared Errors for the Generalized
    %   Voigt Model according to the Lee and Radok indentation
    %   configuration, given a set of input parameters (params).
    %   This function performs fitting for all force curves
    %   separately, by taking in an additional index (idx) compared
    %   to the standard SSE_Voigt function. This is intended
    %   solely for Force Map analysis, wherein each pixel
    %   (containing a single force curve) is treated for analysis.
    %   The output is the same as for the SSE_Voigt function - a
    %   Sum of Squared Errors for that particular pixel.

    % Get Z-Transform Data from Curve
    [Q_hz,F_t,~,h_t,~,~,alphaInit] = zTransformCurve(data,smoothOpt,windowsize,thinPixel);
    U_hz = 1./(Q_hz);
    
    % Get functions
    n_terms = (numel(params)-2)/2;
    [U_storage,U_loss] = getZModels('voigt',n_terms,alphaInit);
    
    % Calculate harmonics
    omega = linspace(-pi,pi,numel(U_hz));
    Uloss_hz = abs(imag(U_hz(omega>=0)));
    Ustorage_hz = real(U_hz(omega>=0));
    omega = omega(omega>=0);
    
    % Trim negative (unphysical) numbers
    inds = find(or(Uloss_hz<0,Ustorage_hz<0));
    omega(inds) = [];
    Uloss_hz(inds) = [];
    Ustorage_hz(inds) = [];
    
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
%     plot(omega,Ustorage_hz,'b-','linewidth',3)
%     hold on
%     plot(omega,U_storage(params,omega),'g-','linewidth',3)
%     xlabel('Angle [radians]')
%     ylabel('Storage Compliance [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     nexttile
%     plot(omega,Uloss_hz,'b-','linewidth',3)
%     hold on
%     plot(omega,U_loss(params,omega),'g-','linewidth',3)
%     xlabel('Angle [radians]')
%     ylabel('Loss Compliance [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     disp('pause');
    
    % calculate global residual
    SError = @(ydata,ymodel) sqrt( sum( ((ydata-ymodel).^2) ./ (numel(ydata) - 2), 'all') ); % Standard Error
%     sse = sum((vertcat(Ustorage_hz,Uloss_hz)-vertcat(U_storage(c,omega),U_loss(c,omega))).^2,'all');
    sse = SError(vertcat(Ustorage_hz,Uloss_hz),vertcat(U_storage(params,omega),U_loss(params,omega)));
    
    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
    
end % End Voigt SSE Z-Transform Map Function