function [] = clusteringAccuracyCalc(originalPath,saveLabel)
%CLUSTERINGACCURACYCALC Calculate Accuracy of Clustering for Simulated Data
%   This function takes in a path containing many subfolders which will be
%   recursively analyzed to determine the accuracy of each clustering
%   method performed. The goal is to summarize the results from many
%   clustering runs on simulated data to determine which observable works
%   best for various noise levels.

% Check to see if there are subdirectories
dirContents = dir(originalPath);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root filesep sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
else
    Folders = {originalPath};
end
    
% Make our output struct
outStruct = struct;

% Begin looping through the directories or files
for i_dir = 1:length(Folders)
        
    % Issue wrapper which allows script to continue if there is an empty
    % directory, or other issue during processing.
    try
    
        fprintf('\nProcessing Directory %d of %d\n\n',i_dir,length(Folders));
        
        path = Folders{i_dir};
        Files = dir([path filesep '*Results*zTransform*.mat']);
        
        for j_dir = 1:length(Files)
            resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');
%             resultsStruct = matfile([Files(j_dir).folder filesep Files(j_dir).name]);
            
            varNames = fields(resultsStruct);

            temp = strsplit(path,filesep);
            dirLabel = temp{end};
            if ~isfield(outStruct, dirLabel)
                % Create a placeholder struct where we store all of the
                % cluster results for each type of analysis.
                temp = struct;
                for ii = 1:6
                    switch ii
                        case 1
                            temp(ii).clusterVar = 'force';

                        case 2
                            temp(ii).clusterVar = 'indentation';

                        case 3
                            temp(ii).clusterVar = 'storage';

                        case 4
                            temp(ii).clusterVar = 'loss';

                        case 5
                            temp(ii).clusterVar = 'angle';

                        case 6
                            temp(ii).clusterVar = 'relaxance';                                
                    end
                    temp(ii).clusterMap2D = [];
                    temp(ii).lastUpdate = '';
                    temp(ii).clusterAccuracy = [];
                end

                outStruct.(dirLabel).clusterData = temp;

            end
            
            for j = 1:numel(varNames)
                
                if ~contains(varNames{j},'zTransform')
                    continue;
                end
                
                if ~isfield(resultsStruct.(varNames{j}), 'trueBinsMap')
                    % This is not a file we want to use for accuracy
                    % calculation because there is no "ground truth". Skip
                    % this variable.
                    continue;
                end
                
                % This should only run once to store the ground truth bins
                if ~isfield(outStruct, 'trueBinsMap')
                    outStruct.trueBinsMap = resultsStruct.(varNames{j}).trueBinsMap;
                end
                
                nBins = max(unique(outStruct.trueBinsMap),[],'all');
                clusterAcc = 0;
                
                for jj = 1:numel([resultsStruct.(varNames{j}).clusterData(:)])

                    if isempty(resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D)
                        % We have not analyzed this type of noise. Continue
                        % to the next type.
                        continue;
                    end
                    
                    binNums = perms(1:nBins);
                    
                    for kk = 1:size(binNums,1)
                        
                        % The clustered bins might not be in the correct
                        % "orientation". This would mean that the clustering
                        % could have correctly separated the regions, BUT it
                        % assigned the innermost bin as #1 instead of #3 for
                        % example. So we will loop through the combos and
                        % ensure we have the highest accuracy calculation
                        % possible.
                        
                        % Cycle through the bin orientations
                        outputMap = resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D;
                        mapDataClusters = NaN(size(outputMap));
                        for kbin = 1:size(binNums,2)
                            mapDataClusters(outputMap == kbin) = binNums(kk,kbin);
                        end
                        
                        binDelta = (mapDataClusters ~= resultsStruct.(varNames{j}).trueBinsMap);
                        temp = 1 - (sum(binDelta,'all') / numel(binDelta));
                        
                        if kk == 1
                            
                            outStruct.(dirLabel).clusterData(jj).clusterMap2D = resultsStruct.(varNames{j}).clusterData(jj).clusterMap2D;
                            outStruct.(dirLabel).clusterData(jj).clusterAccuracy = resultsStruct.(varNames{j}).clusterData(jj).clusterAccuracy;
                            outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);
                            
                        elseif temp > outStruct.(dirLabel).clusterData(jj).clusterAccuracy/100
                            
                            clusterAcc = temp;
                            outStruct.(dirLabel).clusterData(jj).clusterMap2D = mapDataClusters;
                            outStruct.(dirLabel).clusterData(jj).clusterAccuracy = 100*clusterAcc;
                            outStruct.(dirLabel).clusterData(jj).lastUpdate = datestr(now);
%                             fprintf('\nThe Clustering Accuracy (%s) was %.2f%%\n',resultsStruct.(varNames{j}).clusterData(jj).clusterVar,100*clusterAcc);
                        
                        end
                        
                    end
                    
                end
                
                % Now, determine the best observable
                accData = [outStruct.(dirLabel).clusterData(:).clusterAccuracy];
                accData(find(strcmp({outStruct.(dirLabel).clusterData(:).clusterVar},'indentation'),1)) = NaN;
                [~,bestid] = max(accData,[],'omitnan');
                [~,worstid] = min(accData,[],'omitnan');
                
                outStruct.(dirLabel).bestTarget = outStruct.(dirLabel).clusterData(bestid).clusterVar;
                outStruct.(dirLabel).bestAccuracy = outStruct.(dirLabel).clusterData(bestid).clusterAccuracy;
                outStruct.(dirLabel).worstTarget = outStruct.(dirLabel).clusterData(worstid).clusterVar;
                outStruct.(dirLabel).worstAccuracy = outStruct.(dirLabel).clusterData(worstid).clusterAccuracy;
                
                fprintf('\nResults for %s:\n\n',dirLabel);
                fprintf('Best Accuracy: %f%%, %s\n',outStruct.(dirLabel).bestAccuracy,outStruct.(dirLabel).bestTarget);
                fprintf('Worst Accuracy: %f%%, %s\n\n',outStruct.(dirLabel).worstAccuracy,outStruct.(dirLabel).worstTarget);
                
            end
            
        end
        
        clearvars resultsStruct varNames dirLabel
        
    catch ERROR
        
        fprintf('ERROR Clustering Directory #%d of %d\n',i_dir,length(Folders));
        fprintf('The identifier was:\n%s',ERROR.identifier);
        fprintf('Message:%s\n',ERROR.message);
        fprintf('Line Number:%d\n',ERROR.stack(end).line);
        fprintf('Skipping to next directory...\n');
        
    end
    
end

% Save the summary file
save([originalPath filesep 'ClusterSummaryData_' saveLabel],'-struct','outStruct','-v7.3');

end

