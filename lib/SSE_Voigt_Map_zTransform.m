function sse = SSE_Voigt_Map_zTransform(data,params,smoothOpt,windowsize,elasticSetting,fluidSetting)
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
    [Q_hz,~,~,~,~,f_hz,alphaInit] = zTransformCurve(data,smoothOpt,windowsize);
    U_hz = 1./(Q_hz);
    
    % Define Equations
    omega = (2.*pi.*f_hz);

    % Get functions
    n_terms = (numel(params)-2)/2;
    [U_storage,U_loss] = getZModels('voigt',n_terms,alphaInit);
    
    % Calculate harmonics
    Uloss_hz = abs(imag(U_hz));
    Ustorage_hz = real(U_hz);
    
    % calculate global residual
    sse = sum((vertcat(Ustorage_hz,Uloss_hz)-vertcat(U_storage(c,omega),U_loss(c,omega))).^2,'all');

    [tauInds,modulusInds] = getParamIndices(params);
    ub = zeros(size(params))+eps;
    lb = zeros(size(params));

    ub(modulusInds) = 1e2;
    lb(modulusInds) = 1e-12;

    tauCenters = data{8}.*(10.^( (1:length(params(3:2:end)))-1 ));
    ub(tauInds) = tauCenters*10;
    lb(tauInds) = tauCenters/10;

    if length(params) > 1
        if fluidSetting
            ub(2) = 1;
            lb(2) = 0;
        end
    end

    if any(ub-params < 0) || any(params-lb < 0)
        sse = Inf;
    end
    
end % End Voigt SSE Z-Transform Map Function