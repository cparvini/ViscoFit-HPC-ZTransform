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
                case 5
                    thinSample = varargin{i};
                case 6
                    h_finite = varargin{i};
            end
        end
    else
        nu = varargin;
    end
end

c = NaN;
beta = NaN;
switch tipGeom
    case "spherical"
        c = (8*sqrt(tipSize))./(3*(1-nu));
        beta = 1.5;
    case "conical"
        c = (2.*tan(tipSize.*pi./180))./(pi.*(1-nu.^2));
        beta = 2;
    case "punch"
        % Lopez-Guerra & Solares (Beilstein J. of Nanotech., 2017)
        c = (4*tipSize)./(1-nu);
        beta = 1;
    case "pyramidal"
        % Weber, Iturri, Benitez, and Toca-Herrera (Microscopy Research and Technique, 2019)
        % Modified Sneddon
        c = (tan(tipSize.*pi./180))./(sqrt(2).*(1-nu.^2));
        beta = 2;
end

BEC = 1;
if thinSample
    switch tipGeom
        case "spherical"
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            a = sqrt(tipSize.*indentation);
            BEC = (1/(h_finite^0)) + ((1.133.*(a))./(h_finite)) + ...
                ((1.497.*(a.^2))./(h_finite.^2)) + ...
                ((1.469.*(a.^3))./(h_finite.^3)) + ...
                ((0.755.*(a.^4))./(h_finite.^4));

        case "conical"
            % Defined per Garcia & Garcia (Biophy. Journ., 2018)
            a = tan(tipSize.*pi./180);
            BEC = (1/(h_finite^0)) + ((0.721.*(indentation).*(a))./(h_finite)) + ...
                ((0.650.*(indentation.^2).*(a.^2))./(h_finite^2)) + ...
                ((0.491.*(indentation.^3).*(a.^3))./(h_finite^3)) + ...
                ((0.225.*(indentation.^4).*(a.^4))./(h_finite^4));

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
            BEC = (1/(h_finite^0)) + ((0.721.*(indentation).*(a))./(h_finite)) + ...
                ((0.650.*(indentation.^2).*(a.^2))./(h_finite^2)) + ...
                ((0.491.*(indentation.^3).*(a.^3))./(h_finite^3)) + ...
                ((0.225.*(indentation.^4).*(a.^4))./(h_finite^4));

    end
end

% Make our Dirac Delta array for the elastic term
diracArray = zeros(size(time));
diracArray((time-dt)<2*eps) = 1;
    
if length(params) > 2
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

convData = zeros(size(diracArray));
for i = 1:length(startInd)
    temp = convnfft(indentation(startInd(i):endInd(i)).^(beta),Q(startInd(i):endInd(i)),'full');
    convData(startInd(i):endInd(i)) = temp(1:(1+endInd(i)-startInd(i)));
end

out = (c.*convData.*dt)./BEC;

% Release Variables Manually
temp = [];
convData = [];
diracArray = []; 
Q = [];
Q_arms = [];
time_mat = [];

end