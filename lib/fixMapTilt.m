function [heightOut] = fixMapTilt(mapSizeList,pixelHeightList,zeroSubstrate,varargin)
%fixMapTilt Correct for tilted substrate on QI maps
%   Perform a plane-fit on the outer boundary of a QI map. This function
%   will use a polynomial plane-fit to flatten substrates and remove any
%   extraneous features from the boundary (so they don't throw off the tilt
%   of the image). Additional options include the defined fit order
%   (default to 1st-order) and a flag for iterative fit comparisons which
%   will optimize the apparent error of increasing polynomial orders to
%   determine which order polynomial is correct for a given image.

fitOrder = 1;
optimizeOrder = false;
if nargin > 3
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        fitOrder = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        optimizeOrder = varargin{i};                        
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

% Includes substrate pixels
Mask2D(Z - min(Z,[],'all') < 100e-9) = true;

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

% Now, screen the heights for outliers as well
ids = isoutlier(z);
x(ids) = [];
y(ids) = [];
z(ids) = [];

% Initialize Correction Variables
zShift = [];
x2 = reshape(X,numel(pixelHeightList),1);
y2 = reshape(Y,numel(pixelHeightList),1);

if optimizeOrder

    % This will iterate through all of the possible fit cases and return
    % the height correction for each pixel and the polynomial order
    % recommended.
    fitRes = Inf;
    fitOrd = NaN;
    for ii = 1:4
        
        tempRes = Inf;
        loopOrd = num2str(ii);
        [sf,gof] = fit([x, y],z,['poly' loopOrd loopOrd]);
        zShiftTemp = -feval(sf,[x2,y2]);
%         tempRes = gof.rmse;
        tempRes = gof.adjrsquare;
                
        % If this iteration improves the error, save the zShift
        if tempRes < fitRes
            zShift = zShiftTemp;
            fitOrd = ii;
        end
        
    end
    
    fprintf('\nThe optimal substrate-leveling polynomial was of order-%d\n',fitOrd);
    
else
    
    loopOrd = num2str(fitOrder);
    [sf,~] = fit([x, y],z,['poly' loopOrd loopOrd]);
    zShift = -feval(sf,[x2,y2]);
    
end

if zeroSubstrate
    heightOut = cell2mat(pixelHeightList)+zShift';
    heightOut = num2cell(heightOut-min(heightOut,[],'all'));
else
    heightOut = num2cell(cell2mat(pixelHeightList)+zShift');
end

% % Test Plot
% figure
% scatter3(x,y,z,'bo')
% hold on
% scatter3(x(isoutlier(z)),y(isoutlier(z)),z(isoutlier(z)),'rx')
% surf(X,Y,reshape(feval(sf,[reshape(X,numel(pixelHeightList),1),reshape(Y,numel(pixelHeightList),1)]),size(Z)))
% % surf(X,Y,reshape(cell2mat(pixelHeightList)+zShift',size(Z)))
% % surf(X,Y,reshape(cell2mat(heightOut),size(Z)))
% hold off

end

