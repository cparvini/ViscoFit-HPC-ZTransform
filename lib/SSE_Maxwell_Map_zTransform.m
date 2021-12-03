function sse = SSE_Maxwell_Map_zTransform(data,params,ub,lb,smoothOpt,windowsize,thinPixel,elasticSetting,fluidSetting)
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
    [Q_hz,F_t,~,h_t,~,~,alphaInit] = zTransformCurve(data,smoothOpt,windowsize,thinPixel);

    % Get functions
    n_terms = max((numel(params)-2)/2,0);
    [Q_func,Q_storage,Q_loss] = getZModels('maxwell',n_terms,alphaInit);
    
    % Calculate harmonics
    omega = linspace(-pi,pi,numel(Q_hz));
    Qloss_hz = abs(imag(Q_hz(omega>=0)));
    Qstorage_hz = real(Q_hz(omega>=0));
    Q_hz = Q_hz(omega>=0);
    omega = omega(omega>=0);
    
    % Trim negative (unphysical) numbers
    inds = find(or(or(Q_hz<0,Qloss_hz<0),Qstorage_hz<0));
    omega(inds) = [];
    Q_hz(inds) = [];
    Qloss_hz(inds) = [];
    Qstorage_hz(inds) = [];
    
    % Error weighting (if desired)
%     weightsHz = ones(size(omega));
    weightsHz = normpdf(linspace(0,10,numel(omega)),0,0.5);
%     weightsT = ones(size(F_t));
    weightsT = flip(normpdf(linspace(0,5,numel(F_t)),0,1.5));

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
%     plot(omega,Qstorage_hz,'b-','linewidth',3)
%     hold on
%     plot(omega,Q_storage(params,omega),'g-','linewidth',3)
%     xlabel('Angle [radians]')
%     ylabel('Storage Modulus [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     nexttile
%     plot(omega,Qloss_hz,'b-','linewidth',3)
%     hold on
%     plot(omega,Q_loss(params,omega),'g-','linewidth',3)
%     xlabel('Angle [radians]')
%     ylabel('Loss Modulus [Pa]')
%     legend('Data','Fit')
%     hold off
% 
%     disp('pause');
    
    % calculate global residual
    SError = @(ydata,ymodel,weights) sqrt( sum( weights.*(((ydata-ymodel).^2) ./ (numel(ydata) - 2)), 'all') ); % Standard Error
%     sse = sum(weights.*((vertcat(Qstorage_hz,Qloss_hz)-vertcat(Q_storage(params,omega),Q_loss(params,omega))).^2),'all');
%     sse = SError(vertcat(Qstorage_hz,Qloss_hz),vertcat(Q_storage(params,omega),Q_loss(params,omega)));
    sse = abs(SError(Q_hz,Q_func(params,omega),weightsHz)/rms(Q_hz))+...
        abs(SError(Qloss_hz,Q_loss(params,omega),weightsHz)/rms(Qloss_hz))+...
        abs(SError(Qstorage_hz,Q_storage(params,omega),weightsHz)/rms(Qstorage_hz))+...
        abs(SError(F_t,LR_Maxwell(params,data{1},data{2},...
            h_t,data{5},data{6},data{7},elasticSetting,fluidSetting,...
            data{9},data{10}),weightsT)/rms(F_t));
    
    ids = true(size(params));
    if ~elasticSetting
        ids(1) = false;
    end
    if numel(params) > 1
        if ~fluidSetting
            ids(2) = false;
        end
    end
    
    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
    
end % End Maxwell SSE Z-Transform Map Function