clear all
close all
clc

% Make Maps for Each Noise Type
startDir = 'D:\AFM Data\! Simulation';  % Directory with simulated data
mapSize = [128 128];                    % Output Map Size
noiseMag = [30, 40, 50];                % SNR
saveLabel = {'NoiselessSim',...
    'BrownNoiseSim',...
    'AWGNsim',...
    'PinkNoiseSim'};
noiseType = {'none',...
    'brown',...
    'awgn',...
    'pink'};                            % Noise Color

% Create the maps
for i = 1:numel(noiseMag)
    saveDir = ['D:\AFM Data\! Simulation' filesep...
        'SNR-' num2str(noiseMag(i)) filesep];
    for j = 1:numel(noiseType)
        createTestMap(startDir,saveDir,...
            [saveLabel{j} 'SNR' num2str(noiseMag(i))],...
            mapSize,noiseType{j},noiseMag(i));
    end
end