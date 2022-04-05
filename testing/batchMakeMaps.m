clear all
close all
clc

% Make Maps for Each Noise Type
startDir = 'D:\AFM Data\! Simulation';  % Directory with simulated data
mapSize = [128 128];                    % Output Map Size
noiseMag = [25 45; 30 30; 30 40];       % SNR
saveLabel = {'NoiselessSim',...
    'BrownNoiseSim',...
    'AWGNsim',...
    'PinkNoiseSim'};
noiseType = {'none',...
    'brown',...
    'awgn',...
    'pink'};                            % Noise Color

% Create the maps
for i = 1:size(noiseMag,1)
    saveDir = ['D:\AFM Data\! Simulation' filesep...
        'SNR-' strrep(num2str(noiseMag(i,:)),' ','') filesep];
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    for j = 1:numel(noiseType)
        createTestMap(startDir,saveDir,...
            [saveLabel{j} 'SNR' strrep(num2str(noiseMag(i,:)),' ','')],...
            mapSize,noiseType{j},noiseMag(i,:));
    end
end