function sse = SSE_PLR_Map_zTransform(data,params,ub,lb,smoothOpt,windowsize,thinPixel,varargin)
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

    % Get Z-Transform Data from Curve
    [Q_hz,~,~,~,~,~,alphaInit] = zTransformCurve(data,smoothOpt,windowsize,thinPixel);
    
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
    
    modelfun = @(alpha,x) (1/numel(x)).*fftshift(fft(x.*(exp(-alpha.*(1:numel(x))))));
    E_hz = modelfun(alphaInit,E_t);
    
    % Calculate harmonics
    omega = linspace(-pi,pi,numel(Q_hz));
    Qloss_hz = abs(imag(Q_hz(omega>=0)));
    Qstorage_hz = real(Q_hz(omega>=0));
    E_hz = E_hz(omega>=0);
    
    % Trim negative (unphysical) numbers
    inds = find(or(Qloss_hz<0,Qstorage_hz<0));
    Qloss_hz(inds) = [];
    Qstorage_hz(inds) = [];
    E_hz(inds) = [];
    
    % calculate global residual
    SError = @(ydata,ymodel) sqrt( sum( ((ydata-ymodel).^2) ./ (numel(ydata) - 2), 'all') ); % Standard Error
%     sse = sum((vertcat(Qstorage_hz,Qloss_hz)-vertcat(real(E_hz),abs(imag(E_hz)))).^2,'all');
    sse = SError(vertcat(Qstorage_hz,Qloss_hz),vertcat(real(E_hz),abs(imag(E_hz))));

    
    if any(ub(1:length(params))-params < 0) || any(params-lb(1:length(params)) < 0)
        sse = Inf;
    end
    
end % End PLR SSE Z-Transform Map Function