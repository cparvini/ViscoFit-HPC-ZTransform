function [out] = LR_Voigt(params,time,dt,force,tipSize,varargin)
%LR_Voigt Calculate Generalized Kelvin-Voigt L&R Action Integral
%   This function generates the Generalized K-V Model response to the
%   force (F) that is applied by a spherical indenter in time (time) 
%   with tip radius (radius). This model is an application of the Lee and
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

% Calculate coefficient for the action integral
switch tipGeom
    case "spherical"
        c = (3*(1-nu))./(8*sqrt(tipSize));
    case "conical"
        c = 1./((2.*tan(tipSize.*pi./180))./(pi.*(1-nu.^2)));
end

% Make our Dirac Delta array for the elastic term
diracArray = zeros(size(time));
diracArray(time-dt<2*eps) = 1;

if length(params) > 1
    % Make our time matrix (for all the arms)
    time_mat = row2mat(time,length(params(3:2:end)));

    % Calculate the amount of relaxation that occurs for all arms in time. This
    % quantity in time will be removed from the initial modulus response, Eg
    U_arms = sum(params(3:2:end)./params(4:2:end).*exp(-time_mat./params(4:2:end)),1);

    % Add the initial compliance, Jg, such that it is increased to Je as time
    % stretches toward infinity
    if fluidSetting
        % Steady-state fluidity is added to the Retardance
        U = params(1).*(diracArray./dt) + U_arms + params(2);
    else
        U = params(1).*(diracArray./dt) + U_arms;
    end
else
    U = params(1).*(diracArray./dt);
end

% Calculate the action integral quantity
startInd = find(diracArray);
endInd = horzcat(find(diracArray)-1,numel(diracArray));
endInd(1) = [];

convData = zeros(size(diracArray));
for i = 1:length(startInd)
    temp = convnfft(force(startInd(i):endInd(i)),U(startInd(i):endInd(i)),'full');
    convData(startInd(i):endInd(i)) = temp(1:(1+endInd(i)-startInd(i)));
end

out = c.*convData.*dt;

end