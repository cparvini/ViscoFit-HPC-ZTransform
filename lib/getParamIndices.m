function [tauInds,modulusInds] = getParamIndices(beta)
%getParamIndices Get the Indices for All Moduli and Characteristic Times
%   This function creates a list of the indices for the Moduli and
%   Characteristic times inside an array of size beta.

% Ignore the SS fluidity, which will not have the same bounds as normal 
% characteristic times
tauInds = (4:2:length(beta));

% Incude all moduli, because even the elastic/glassy stiffness must fall
% within reasonable bounds.
modulusInds = (1:2:length(beta));

end

