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
percentImprove = 0.1;
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
                case 3
                    if ~isempty(varargin{i})
                        percentImprove = varargin{i}/100;                        
                    end
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

mapSize = mapSizeList{1};

if iscell(mapSize)
    while iscell(mapSize)
        mapSize = mapSize{1};
    end
end

% axes meshgrid for scattering data
em1 = floor(mapSize(1)/20);
em2 = floor(mapSize(2)/20);
% em1 = 1;
% em2 = 1;
xdata = 1:mapSize(1);
ydata = flip(1:mapSize(2));
[X, Y] = meshgrid(xdata,ydata);
Z = reshape([pixelHeightList{:}],mapSize);
Mask2D = false(size(Z));
Mask2D(:,[1:em2,(end-em2):end]) = true;
Mask2D([1:em1,(end-em1):end],:) = true;

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
ids = isoutlier(z,'percentiles',[5 95]);
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
    fitOrd = {NaN,NaN};
    for ii = 1:4
        
        tempRes = Inf;
        
        for jj = 1:4
            [sf,gof] = fit([x, y],z,['poly' num2str(ii) num2str(jj)]);
            zShiftTemp = -feval(sf,[x2,y2]);
    %         tempRes = gof.rmse;
            tempRes = gof.adjrsquare;

            % If this iteration improves the error, save the zShift
            if tempRes < fitRes*(1-percentImprove)
                zShift = zShiftTemp;
                fitOrd = {ii,jj};
            end
        end
        
    end
    
    fprintf('\nThe optimal substrate-leveling plane is...\na %d-order polynomial in X, and\na %d-order polynomial in Y\n\n',fitOrd{1},fitOrd{2});
    
else
    
    if any(fitOrder <= 0) || (numel(fitOrder) < 1)
        fitOrder = 1;
    end
    
    if numel(fitOrder) == 1
        loopOrd = {num2str(fitOrder),num2str(fitOrder)};
    elseif numel(fitOrder) == 2
        loopOrd = {num2str(fitOrder(1)),num2str(fitOrder(2))};
    else
        loopOrd = {num2str(fitOrder(1)),num2str(fitOrder(2))};
    end
    
    [sf,~] = fit([x, y],z,['poly' loopOrd{1} loopOrd{2}]);
    zShift = -feval(sf,[x2,y2]);
    
end

if zeroSubstrate
    heightOut = cell2mat(pixelHeightList)+zShift';
%     heightOut = num2cell(heightOut-min(heightOut,[],'all'));
    heightOut = num2cell(heightOut-min(heightOut(heightOut>0),[],'all'));
else
    heightOut = num2cell(cell2mat(pixelHeightList)+zShift');
end

% % Test Plot
% figure
% % scatter3(x,y,z,'bo')
% hold on
% % scatter3(x(isoutlier(z)),y(isoutlier(z)),z(isoutlier(z)),'rx')
% surf(X,Y,reshape(feval(sf,[reshape(X,numel(pixelHeightList),1),reshape(Y,numel(pixelHeightList),1)]),size(Z)))
% % surf(X,Y,reshape(cell2mat(pixelHeightList)+zShift',size(Z)))
% % surf(X,Y,reshape(cell2mat(heightOut),size(Z)))
% hold off

end