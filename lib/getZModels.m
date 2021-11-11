function [storageFunc,lossFunc] = getZModels(model,n_terms,alphaInit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch lower(model)
    case 'maxwell'
        if n_terms > 0
%             storageFunc = @(c,omega) sum(c(1:2:end)) - ...
%                 sum( (c(3:2:end).*repmat(ones(size(omega)),n_terms,1).*(((1+(c(4:2:end)./dt).*alphaInit).^2) ./ ...
%                 ((c(4:2:end)./dt).^2))) ... 
%                 ./ ( ((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)) .* ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1)))...
%                 , 1);
            storageFunc = @(c,omega) sum(c(1:2:end)) - ...
                sum( ( (c(3:2:end).*repmat(ones(size(omega)),n_terms,1)) ./ ((1 + alphaInit.*((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1)) ) ... 
                .* ( ((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                ./ (((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                + repmat(omega.^2,n_terms,1)) )...
                , 1);
%             lossFunc = @(c,omega) sum( repmat(omega,n_terms,1).*(c(3:2:end).*(c(4:2:end)./dt).*repmat(ones(size(omega)),n_terms,1) .* ...
%                 (((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2))) ./ ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2).*repmat(ones(size(omega)),n_terms,1)) .* ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1)))...
%                 , 1);
            lossFunc = @(c,omega) sum(c(1:2:end)) - ...
                sum( repmat(omega,n_terms,1) ...
                .* ( (c(3:2:end).*c(4:2:end).*repmat(ones(size(omega)),n_terms,1)) ./ ((1 + alphaInit.*((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1)) ) ... 
                .* ( ((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                ./ (((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                + repmat(omega.^2,n_terms,1)) )...
                , 1);
        else
            storageFunc = @(c,omega) c(1).*ones(size(omega));
            lossFunc = @(c,omega) zeros(size(omega));
        end
        
    case 'voigt'
        if n_terms > 0
%             storageFunc = @(c,omega) sum(c(1:2:end)) - ...
%                 sum( (c(3:2:end).*repmat(ones(size(omega)),n_terms,1).*(((1+(c(4:2:end)./dt).*alphaInit).^2) ./ ...
%                 ((c(4:2:end)./dt).^2))) ... 
%                 ./ ( ((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)) .* ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1)))...
%                 , 1);
            storageFunc = @(c,omega) c(1) - ...
                sum( ( (c(3:2:end).*repmat(ones(size(omega)),n_terms,1)) ./ ((1 + alphaInit.*((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1)) ) ... 
                .* ( ((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                ./ (((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                + repmat(omega.^2,n_terms,1)) )...
                , 1);
%             lossFunc = @(c,omega) sum( repmat(omega,n_terms,1).*(c(3:2:end).*(c(4:2:end)./dt).*repmat(ones(size(omega)),n_terms,1) .* ...
%                 (((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2))) ./ ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2).*repmat(ones(size(omega)),n_terms,1)) .* ...
%                 ((((1+(c(4:2:end)./dt).*alphaInit).^2)./((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1) + repmat(omega.^2,n_terms,1)))...
%                 , 1);
            lossFunc = @(c,omega) sum(c(1:2:end)) - ...
                sum( repmat(omega,n_terms,1) ...
                .* ( (c(3:2:end).*c(4:2:end).*repmat(ones(size(omega)),n_terms,1)) ./ ((1 + alphaInit.*((c(4:2:end)./dt).^2)).*repmat(ones(size(omega)),n_terms,1)) ) ... 
                .* ( ((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                ./ (((((1+(c(4:2:end)./dt).*alphaInit).*repmat(ones(size(omega)),n_terms,1)).^2) ./ ((c(4:2:end)./dt).^2).*repmat(ones(size(omega)),n_terms,1)) ...
                + repmat(omega.^2,n_terms,1)) )...
                , 1);
        else
            storageFunc = @(c,omega) c(1).*ones(size(omega));
            lossFunc = @(c,omega) zeros(size(omega));
        end
        
end

end

