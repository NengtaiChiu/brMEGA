function [artifactResidueIndex, spectralConcentration, atM ] = nlem_performance( cleanData, rawData, stimLoc )
% Function to compute artifact residual index and spectral concentration
% Artifact Residual Index
% extract interstim intervals
sampling_rate = 250;
interStimData = cell(length(stimLoc), 2);
interStimData{1, 1} = 1;
for i = 2:length(stimLoc)
    interStimData{i, 1} = i;
    k = stimLoc(i) - stimLoc(i - 1);
    if 0.75 * sampling_rate < k < 1.5 * sampling_rate % 1.5*sampling_rate
        interStimData{i, 2} = rawData(stimLoc(i - 1) + 0.5*sampling_rate:stimLoc(i) -0.2*sampling_rate);
    else
        interStimData{i, 2} = rawData(stimLoc(i) -0.5*sampling_rate :stimLoc(i) -0.2*sampling_rate);
    end
end
interStimData = interStimData(~cellfun(@isempty, interStimData(:, 2)), :);

% find nn closest intervals to each artifact
neighborISData = cell(length(stimLoc), 1);
nn = 3;%
for i = 1:length(stimLoc)
    t = find(cell2mat(interStimData(:, 1)) >= ...
        min(length(stimLoc) - (2 * nn + 1), max(1, i - nn)), 2 * nn + 1, 'first');
    neighborISData{i} = cell2mat(interStimData(t, 2)');
end
% take the average l2 and linfty energies
interStimEnergy2 = zeros(length(stimLoc), 1);
interStimEnergyInf = zeros(length(stimLoc), 1);
for i = 1:length(stimLoc)
    interStimEnergy2(i) = median(abs(neighborISData{i} - median(neighborISData{i})));
    interStimEnergyInf(i) = prctile(abs(neighborISData{i} - median(neighborISData{i})), 95);
end

artifactEnergy2 = zeros(length(stimLoc), 1);
artifactEnergyInf = zeros(length(stimLoc), 1);
r_len = 0.5*sampling_rate;
l_len = 0.2*sampling_rate;
for i = 2:length(stimLoc)-1
    artifactEnergy2(i) = median(abs(cleanData(stimLoc(i) - l_len:stimLoc(i) + r_len) - ...
        median(cleanData(stimLoc(i) - l_len:stimLoc(i) + r_len))));
    artifactEnergyInf(i) = max(abs(cleanData(stimLoc(i) - l_len:stimLoc(i) + r_len) - ...
        median(cleanData(stimLoc(i) - l_len:stimLoc(i) + r_len))));
end

artifactResidueIndex = zeros(length(stimLoc), 1);
for i = 1:length(stimLoc)
    term1 = (1/ 2) * ((artifactEnergy2(i) / interStimEnergy2(i)) + (interStimEnergy2(i) / artifactEnergy2(i)));
    term2 = (1/ 2) * ((artifactEnergyInf(i) / interStimEnergyInf(i)) + (interStimEnergyInf(i) / artifactEnergyInf(i)));
    at(1,i) = term1;
    at(2,i) = term2;
    artifactResidueIndex(i) = term1 * term2;
end
artifactResidueIndex = median(artifactResidueIndex);
atM(1) = median(at(1,:));
atM(2) = median(at(2,:));
% spectral concentration
[powerClean,f] = pwelch(cleanData, 2500,[], 2500, 250);
powerRaw = pwelch(rawData, 2500, [], 2500, 250);

f = f(f>0 & f<=200 );

powerClean = powerClean(f>=1 & f<=200);
powerRaw = powerRaw(f>=1 & f<=200);

stimFreq = median(250./diff(stimLoc));
stimFreqHarmonics = stimFreq*(1:10);
fIndices = [];
for iF = 1:length(stimFreqHarmonics)
    fIndices = cat(1,fIndices,find(f>stimFreqHarmonics(iF)-0.1 & ...
        f<stimFreqHarmonics(iF)+0.1));
end
%keyboard
pChange = powerClean./powerRaw;
spectralConcentration = sum(pChange(fIndices)) / sum(pChange(setdiff(1:length(pChange),fIndices)));%keyboard