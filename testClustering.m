function [] = testClustering(startDir,N_workers,varargin)
%TESTCLUSTERING Loop through clustering options for a directory to see what
%the single-map results look like for each observable
%   Taking in a directory (startDir), this function runs successive
%   clustering attempts on the files in that directory to show what
%   clustering results look like in each of the available cases. Plots are
%   created for each case that can be observed later. The number of
%   parallel workers (N_workers) and frequency to plot results at (evalPt)
%   are passed to the singleMapClustering() function.

% User-Defined Settings
correctTilt = true;
hideSubstrate = true;
zeroSubstrate = true;
optimizeFlattening = false;
fillPixels = true;
logSteps = true;
plotIndentation = true;
evalPt = 1000;
n_reps = 10; % number of clustering replicates
maxK = 10; % Max number of cluster bins
if nargin > 1
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            switch i
                case 1
                    if ~isempty(varargin{i})
                        correctTilt = varargin{i};                        
                    end
                case 2
                    if ~isempty(varargin{i})
                        hideSubstrate = varargin{i};
                    end
                case 3
                    if ~isempty(varargin{i})
                        zeroSubstrate = varargin{i};                        
                    end
                case 4
                    if ~isempty(varargin{i})
                        optimizeFlattening = varargin{i};                        
                    end
                case 5
                    if ~isempty(varargin{i})
                        fillPixels = varargin{i};                        
                    end
                case 6
                    if ~isempty(varargin{i})
                        logSteps = varargin{i};                        
                    end
                case 7
                    if ~isempty(varargin{i})
                        plotIndentation = varargin{i};                        
                    end
                case 8
                    if ~isempty(varargin{i})
                        evalPt = varargin{i};                        
                    end
                case 9
                    if ~isempty(varargin{i})
                        n_reps = varargin{i};                        
                    end
                case 10
                    if ~isempty(varargin{i})
                        maxK = varargin{i};                        
                    end
                otherwise
                    fprintf('Passed additional parameters to fit_map() which were not used.');
            end
        end
    end
end

for ii = 1:6
    
    switch ii
        case 1
            clusterTarget = 'force';
        case 2
            clusterTarget = 'indentation';
        case 3
            clusterTarget = 'storage';
        case 4
            clusterTarget = 'loss';
        case 5
            clusterTarget = 'angle';
        case 6
            clusterTarget = 'relaxance';
    end
    
    fprintf('Performing %s clustering (Method %d of 6)...\n',clusterTarget,ii);
   
    singleMapClustering(startDir,...
        N_workers,...
        clusterTarget,...
        correctTilt,...                 % Correct Tilt (fit a polynomial plane to the height data and fix tilt/bend)
        hideSubstrate,...               % Hide Substrate (remove it from analysis)
        zeroSubstrate,...               % Zero Substrate (shift down to zero)
        optimizeFlattening,...          % Optimize Flattening (test different orders of flattening functions)
        fillPixels,...                  % Fill Pixels (use neighbors to fill bad pixels)
        logSteps,...                    % Log-Scaled Steps (change step size for each order of 10 in frequency)
        plotIndentation,...             % Plot Indentation (Include the indentation map in the output plot)
        evalPt,...                      % Evaluation Frequency (frequency, in Hz, to evaluate functions at for plotting in the figure)
        n_reps,...                      % Number of Replicates (number of times to repeat clustering; higher number -> more confidence)
        maxK)                           % Maximum Number of Bins (for clustering, maximum number of "bins" to place data in; all numbers below this will be tested)
    
    fprintf('%s clustering complete.\n\n',clusterTarget);
    
end

end

