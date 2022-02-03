function [dist] = dtwf(x,y)
%DTWF Wrapper for Dynamic Time Warping Comparison of time-series vectors
%   This function takes in an array x (n by 1) for a single time-series
%   measurement and compares it to all of the time series measurements
%   stored in the cell array y (n by m). The one caveat to this program is
%   the in-built nan-removal. This will allow comparison of time series
%   with different lengths, which are usually NaN-padded during processing.

% Initialize
m = size(y,1);
dist = zeros(m,1);
xfin = x(~isnan(x));

for i = 1:m
    yfin = y(i,:);
    yfin = yfin(~isnan(yfin));
    dist(i) = dtw(xfin,yfin);
end

end