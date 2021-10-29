function out = row2mat(x,n)
%row2mat Vertically Stack Row Vector Many Times
%   This function replicates the array input (x) exactly n times in the 1st
%   dimension (i.e. stacking rows vertically). This function is not truly
%   necessary in Matlab, because of the repmat() functionality built-in,
%   but for code-comparability to Python and Julia it has been included.
out = repmat(x,[n 1]);
end

