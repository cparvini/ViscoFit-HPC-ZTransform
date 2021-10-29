function [out] = LR_Maxwell(params,time,dt,indentation,tipSize,varargin)
%LR_Maxwell Calculate Generalized Maxwell L&R Action Integral
%   This function generates the Generalized Maxwell Model response to an
%   applied spherical indentation (indentation) in time (time) with a probe 
%   of tip radius (radius). This model is an application of the Lee and
%   Radok (1960) indentation configuration, which is well-studied in the
%   viscoelastic parameter extraction literature; therein, a spherical
%   probe is assumed to indent a viscoelastic half-space that is previously
%   undisturbed.

% Check the varargin
nu = 0.5; % Default poisson's ratio of the sample to incompressible (0.5)
tipGeom = "spherical";
elasticSetting = 1;
fluidSetting = 0;
thinSample = false;
h_finite = NaN;
if ~isempty(varargin)
    if length(varargin) > 1
        for i = 1:length(varargin)
            switch(i)
                case 1
                    nu = varargin{i};
                case 2
                    tipGeom = varargin{i};
                case 3
                    elasticSetting = varargin{i};
                case 4
                    fluidSetting = varargin{i};
                case 6
                    thinSample = varargin{i};
                case 7
                    h_finite = varargin{i};
            end
        end
    else
        nu = varargin;
    end
end

if thinSample && isnan(h_finite)
    error('You attempted to enforce finite sample thickness, but did not define the sample thickness a-priori. Ensure that you are passing this value to LR_Maxwell()');
end

% Calculate coefficient for the action integral
switch tipGeom
    case "spherical"
        if ~thinSample
            c = (8*sqrt(tipSize))./(3*(1-nu));
            beta = 1.5;
        else
            % Defined per Garcia & Garcia (Nanoscale, 2018), Table 3
            cTaylor = 1./[(8*(mode(tipSize)^(1/2)))./(3*(1-nu))...
                1.133*(8*(mode(tipSize)^(1)))./(3*(1-nu))./(h_finite)...
                1.497*(8*(mode(tipSize)^(3/2)))./(3*(1-nu))./(h_finite.^2)...
                1.469*(8*(mode(tipSize)^(2)))./(3*(1-nu))./(h_finite.^3)...
                0.755*(8*(mode(tipSize)^(5/2)))./(3*(1-nu))./(h_finite.^4)];
            betaTaylor = [3/2 2 5/2 3 7/2];
        end
    case "conical"
        if ~thinSample
            c = (2.*tan(tipSize.*pi./180))./(pi.*(1-nu.^2));
            beta = 2;
        else
            % Defined per Garcia, Guerrero, & Garcia (Nanoscale, 2020),
            % Supplemental Information
            cTaylor = 1./[8*tan(mode(tipSize).*pi./180)./(3*pi)...
                0.721*8*(tan(mode(tipSize).*pi./180).^2)./(3*(h_finite)*pi)...
                0.650*8*(tan(mode(tipSize).*pi./180).^3)./(3*(h_finite.^2)*pi)...
                0.491*8*(tan(mode(tipSize).*pi./180).^4)./(3*(h_finite.^3)*pi)...
                0.225*8*(tan(mode(tipSize).*pi./180).^5)./(3*(h_finite.^4)*pi)];
            betaTaylor = [2 3 4 5 6];
        end
end

% Make our Dirac Delta array for the elastic term
diracArray = zeros(size(time));
diracArray((time-dt)<2*eps) = 1;
    
if length(params) > 1
    % Make our time matrix (for all the arms)
    time_mat = row2mat(time,length(params(3:2:end)));

    % Calculate the amount of relaxation that occurs for all arms in time. This
    % quantity in time will be removed from the initial modulus response, Eg
    Q_arms = sum(params(3:2:end)./params(4:2:end).*exp(-time_mat./params(4:2:end)),1);

    % Add the initial modulus, Eg, such that it is relaxed to Ee as time
    % stretches toward infinity
    if fluidSetting
        % Eg relaxes in time due to the steady-state viscosity stored inside
        % params(2)
        Q = sum(params(1:2:end)).*(diracArray./dt) - (params(1)./params(2).*exp(-time./params(2))) - Q_arms;
    else
        Q = sum(params(1:2:end)).*(diracArray./dt) - Q_arms; % Relax Eg in time
    end
else
    Q = params(1).*(diracArray./dt);
end

% Calculate the action integral quantity
startInd = find(diracArray);
endInd = horzcat(find(diracArray)-1,numel(diracArray));
endInd(1) = [];
if ~thinSample
    convData = zeros(size(diracArray));
    for i = 1:length(startInd)
        temp = convnfft(indentation(startInd(i):endInd(i)).^(beta),Q(startInd(i):endInd(i)),'full');
        convData(startInd(i):endInd(i)) = temp(1:(1+endInd(i)-startInd(i)));
    end

    out = c.*convData.*dt;    
else
    % Introduce the higher order terms
    out = zeros(size(diracArray));
    for i = 1:numel(cTaylor)
        convData = zeros(size(diracArray));
        for j = 1:length(startInd)
            temp = convnfft(indentation(startInd(j):endInd(j)).^(betaTaylor(i)),Q(startInd(j):endInd(j)),'full');
            convData(startInd(j):endInd(j)) = temp(1:(1+endInd(j)-startInd(j)));
        end

        out = out + cTaylor(i).*convData.*dt;
    end
end

% Release Variables Manually
temp = [];
convData = [];
diracArray = []; 
Q = [];
Q_arms = [];
time_mat = [];
cTaylor = [];
betaTaylor = [];

end