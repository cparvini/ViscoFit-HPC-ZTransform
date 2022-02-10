function [heightOut] = fixMapTilt(mapSizeList,pixelHeightList,zeroSubstrate,varargin)
%fixMapTilt Correct for tilted substrate on QI maps
%   Perform a plane-fit on the outer boundary of a QI map. This function
%   will loop until it has successfully flattened a map and removed any
%   extraneous features from the boundary (so they don't throw off the tilt
%   of the image).

fitOrder = 2;
if nargin > 3
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        fitOrder = varargin{i};                        
                    end
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

mapSize = mapSizeList{1};

% axes meshgrid for scattering data
xdata = 1:mapSize(1);
ydata = flip(1:mapSize(2));
[X, Y] = meshgrid(xdata,ydata);
Z = reshape([pixelHeightList{:}],mapSize);
Mask2D = false(size(Z));
Mask2D(:,[1,end]) = true;
Mask2D([1,end],:) = true;

x = reshape(X,numel(pixelHeightList),1);
y = reshape(Y,numel(pixelHeightList),1);
z = reshape(Z,numel(pixelHeightList),1);
mask2D = reshape(Mask2D,numel(pixelHeightList),1);

x(~mask2D) = [];
y(~mask2D) = [];
z(~mask2D) = [];

% Remove outliers
temp = imgradient(Z,'sobel');
temp(~mask2D) = [];
ids = isoutlier(temp);
x(ids) = [];
y(ids) = [];
z(ids) = [];

switch fitOrder
    case 1
        DM = [x,y,ones(size(z))];
        B = DM\z;
        ZM = B(1)*X+B(2)*Y+B(3)*ones(size(X));

        Z1 = ZM - min(ZM,[],'all');
        ZShift = -Z1;
        zShift = reshape(ZShift,numel(pixelHeightList),1);
    
    case 2
        sf = fit([x, y],z,'poly22');
        x2 = reshape(X,numel(pixelHeightList),1);
        y2 = reshape(Y,numel(pixelHeightList),1);
        zShift = -feval(sf,[x2,y2]);
        
    case 3
        sf = fit([x, y],z,'poly33');
        x2 = reshape(X,numel(pixelHeightList),1);
        y2 = reshape(Y,numel(pixelHeightList),1);
        zShift = -feval(sf,[x2,y2]);
        
    otherwise
        warning('fixMapTilt() correction does not allow fitting surface polynomials greater than 3rd order. Defaulting to 3rd order polynomial in x and y.')
end

if zeroSubstrate
    heightOut = cell2mat(pixelHeightList)+zShift';
    heightOut = num2cell(heightOut-min(heightOut,[],'all'));
else
    heightOut = num2cell(cell2mat(pixelHeightList)+zShift');
end

% % Test performance
% figure
% plot3(x,y,z-min(z,[],'all'),'.')
% hold on
% plot(sf,[x,y],z)
% % meshc(X, Y, ZM)
% % meshc(X, Y, ZShift)
% meshc(X, Y, reshape([heightOut{:}],mapSize))
% hold off
% grid on
% xlabel('x'); ylabel('y'); zlabel('h');
% grid on

end

