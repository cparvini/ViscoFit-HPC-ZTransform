function sse = SSE_PLR_Map_zTransform(data,params,smoothOpt,windowsize,varargin)
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
    [Q_hz,~,~,~,~,~,alphaInit] = zTransformCurve(data,smoothOpt,windowsize);

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
    
    modelfun = @(alpha,x) (1/numel(x)).*fftshift(fft(x.*(exp(-alpha.*(1:numel(x))))));
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