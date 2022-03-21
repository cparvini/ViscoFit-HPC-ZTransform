function [SNRforce_record,SNRind_record] = estimateSNR(originalPath)
%ESTIMATESNR Estimate the SNR for the Force and Indentation of a QI map

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

% Issue wrapper which allows script to continue if there is an empty
% directory, or other issue during processing.
try

    % Begin looping through the directories or files
    for i_dir = 1:length(Folders)
        path = Folders{i_dir};
        Files = dir([path filesep '*Results*zTransform*.mat']);

        if isempty(Files)
            error('The directory you selected does not contain a Z-Transform QI map. Please verify your FitResults file is in that directory and the filename contains "zTransform".');
        end

        for j_dir = 1:length(Files)
            resultsStruct = load([Files(j_dir).folder filesep Files(j_dir).name],'-mat');

            varNames = fields(resultsStruct);

            for j = 1:numel(varNames)
                
                if ~contains(varNames{j},'zTransform')
                    continue;
                end
        
                SNRforce = NaN(numel(resultsStruct.(varNames{j}).ViscoClass.forces_cell),1);
                SNRind = NaN(numel(resultsStruct.(varNames{j}).ViscoClass.indentations_cell),1);
                
                for i = 1:numel(resultsStruct.(varNames{j}).ViscoClass.forces_cell)
                    
                    % Quick Linear Fit
                    ydata = resultsStruct.(varNames{j}).ViscoClass.forces_cell{i};
                    xdata = resultsStruct.(varNames{j}).ViscoClass.times_cell{i};
                    
                    if isempty(xdata) || isempty(ydata)
                        continue;
                    end
                    
                    if isrow(ydata) ydata = ydata'; end
                    if isrow(xdata) xdata = xdata'; end
                    xdata = xdata - xdata(1);
                    ydata = ydata - ydata(1);

                %     b = [ones(size(xdata)) xdata] \ ydata;
                %     model = b(1) + b(2)*xdata;

                %     b = xdata \ ydata;
                %     model = b*xdata;

                    p = polyfit(xdata,ydata,2);
                    model = polyval(p,xdata);

%                     % plot
%                     figure(1)
%                     plot(xdata,ydata,'r-',xdata,model,'b--')

                    residual_noise = ydata - model; 
                    snr_calc = mean(ydata.^2) / mean(residual_noise.^2); 
                    SNRforce(i) = 10 * log10(snr_calc);
                %     SNRforce(i) = snr(ydata,model);

                    % Quick Linear Fit
                    ydata = resultsStruct.(varNames{j}).ViscoClass.indentations_cell{i};
                    if isrow(ydata) ydata = ydata'; end
                    ydata = ydata - ydata(1);

                %     b = [ones(size(xdata)) xdata] \ ydata;
                %     model = b(1) + b(2)*xdata;

                %     b = xdata \ ydata;
                %     model = b*xdata;

                    p = polyfit(xdata,ydata,2);
                    model = polyval(p,xdata);

%                     % plot
%                     figure(2)
%                     plot(xdata,ydata,'r-',xdata,model,'b--')

                    residual_noise = ydata - model; 
                    snr_calc = mean(ydata.^2) / mean(residual_noise.^2); 
                    SNRind(i) = 10 * log10(snr_calc);
                %     SNRind(i) = snr(ydata,model);

                end

                fprintf('Resuls for %s...\n',Files(j_dir).name)
                fprintf('The average Force SNR was: %3.2f\n',mean(SNRforce,'omitnan'))
                fprintf('The average Indentation SNR was: %3.2f\n\n',mean(SNRind,'omitnan'))

            end
            
            if ~exist('SNRforce_record','var')
                SNRforce_record = SNRforce;
                SNRind_record = SNRind;
            else
                SNRforce_record = vertcat(SNRforce_record,SNRforce);
                SNRind_record = vertcat(SNRind_record,SNRind);
            end
            
        end
        
    end

catch ERROR
        
    fprintf('ERROR Clustering Directory #%d of %d\n',i_dir,length(Folders));
    fprintf('The identifier was:\n%s',ERROR.identifier);
    fprintf('Message:%s\n',ERROR.message);
    fprintf('Line Number:%d\n',ERROR.stack(end).line);
    fprintf('Skipping to next directory...\n');
    
end

end