function [heightOut] = fixMapTilt(mapSizeList,pixelHeightList,zeroSubstrate)
%fixMapTilt Correct for tilted substrate on QI maps
%   Detailed explanation goes here

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

DM = [x,y,ones(size(z))];
B = DM\z;
ZM = B(1)*X+B(2)*Y+B(3)*ones(size(X));

Z1 = ZM - min(ZM,[],'all');
ZShift = -Z1;
zShift = reshape(ZShift,numel(pixelHeightList),1);

if zeroSubstrate
    heightOut = cell2mat(pixelHeightList)+zShift';
    heightOut = num2cell(heightOut-min(heightOut,[],'all'));
else
    heightOut = num2cell(cell2mat(pixelHeightList)+zShift');
end

% % Test out the linear regression
% figure
% plot3(x,y,z,'.')
% hold on
% meshc(X, Y, ZM)
% meshc(X, Y, ZShift)
% meshc(X, Y, reshape([heightOut{:}],mapSize))
% hold off
% grid on
% xlabel('x'); ylabel('y'); zlabel('h');
% grid on

end

