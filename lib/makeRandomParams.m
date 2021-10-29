function [beta0] = makeRandomParams(oldParams,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,varargin)
%makeRandomParams Create Random Parameters Within Bounds
%   This function takes in the relevant elastic and fluid settings for the
%   viscoelastic model, including the limits within which random guesses
%   should be generated, and provides a random starting point for the
%   fitting to take. By passing the additional argument "newInds", this
%   will specifically replace only the NEW elements inside of the parameter
%   set. This feature is used exclusively to iteratively introduce
%   terms---by feeding forward the old results in "oldParams" and only
%   giving a new starting position for the indices that have never been fit
%   before.

    % Initialize
    newInds = (false(size(ub_randLimit)));
    
    % Checking varargin structure
    if ~isempty(varargin)
        if length(varargin) > 1
            for i = 1:length(varargin)
                switch i
                    case 1
                        newInds = (cell2mat(varargin{i}));
%                     case 2
%                         % new varargin argument....
                    otherwise
                       error('Too many arguments passed to makeRandomParams via varargin!')
                end
            end
        else
            newInds = (cell2mat(varargin));
        end
    end
    
    beta0 = oldParams;

    beta0_temp = zeros(size(ub_randLimit));
    beta0_temp([1 2]) = 10.^(((ub_randLimit([1 2])-lb_randLimit([1 2])).*rand(1,2)+lb_randLimit([1 2])));
    beta0_temp(3:2:end) = 10.^(((ub_randLimit(3:2:end)-lb_randLimit(3:2:end)).*rand(size(ub_randLimit(3:2:end))) + lb_randLimit(3:2:end)));
    beta0_temp(4:2:end) = 10.^((ub_randLimit(4:2:end)-lb_randLimit(4:2:end)).*rand(size(ub_randLimit(4:2:end))) + lb_randLimit(4:2:end));
    
    if ~elasticSetting
        beta0_temp(1) = 0;
    end
    
    if ~fluidSetting
        beta0_temp(2) = 0;
    end
    
    if any(newInds)
        beta0(newInds) = beta0_temp(newInds);
    else
        beta0 = beta0_temp;
    end
    
end

