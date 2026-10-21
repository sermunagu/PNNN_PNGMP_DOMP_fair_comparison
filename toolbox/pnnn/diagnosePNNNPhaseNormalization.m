function report = diagnosePNNNPhaseNormalization(modelFile, measurementFile, degrees)
% Read-only model diagnosis. The zero-channel intervention is NOT a trained model.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
saved = load(modelFile, 'bundle');
bundle = saved.bundle;
capture = load(measurementFile, 'x', 'y');
u = double(capture.x(:));
y = double(capture.y(:));
if bundle.config.pnnn.removeDC, y = y - mean(y); end
model = bundle.models{3};
f = model.featureConfig;
stats = model.normalization;
channel = f.M + 2;
assert(string(f.featMode) == "full", 'This diagnosis addresses the full-mode current imaginary tap.');
id = bundle.split.identificationIndices;
evaluation = setdiff((1:numel(u)).', id);
if nargin < 3, degrees = [0 30 60 90 135 180 270]; end
assert(degrees(1) == 0, 'The reference rotation must be zero.');
learnables = model.network.Learnables;
weightRow = find(learnables.Layer == "fc1" & learnables.Parameter == "Weights");
weights = double(extractdata(learnables.Value{weightRow}));
active = logical(model.masks{weightRow}(:,channel));
rotationTable = table();
featureStats = table();
for index = 1:numel(degrees)
    phase = exp(1j*degrees(index)*pi/180);
    rotatedInput = phase*u;
    predicted = applyDPDModels(bundle, rotatedInput);
    x = rotatedInput;
    if bundle.config.pnnn.removeDC, x = x - mean(x); end
    [features, ~, rotation] = buildPhaseNormDataset(x, zeros(size(x)), f.M, f.orders, f.featMode, f);
    normalized = (features - stats.muX.') ./ stats.sigmaX.';
    rawChannel = features(channel,:).';
    normalizedChannel = normalized(channel,:).';
    if index == 1
        baseline = predicted; % The required first angle is zero; reuse its inference.
        referenceChannel = normalizedChannel;
        featureScale = table((1:size(features,1)).', mean(features,2), std(features,0,2), ...
            max(abs(features),[],2), stats.muX.', stats.sigmaX.', ...
            sqrt(mean(normalized.^2,2)), max(abs(normalized),[],2), ...
            'VariableNames', {'Feature','Mean','Std','MaxAbs','StoredMean','StoredSigma','NormalizedRMS','NormalizedMaxAbs'});
        featureScale.SigmaOverInputRMS = featureScale.StoredSigma/sqrt(mean(abs(x).^2));
        featureScale.StdOverFeatureRMS = featureScale.Std ./ max(sqrt(mean(features.^2,2)), realmin);
        contribution = normalizedChannel * weights(:,channel).';
        other = normalized.' * weights.' - contribution;
        connections = table((1:size(weights,1)).', weights(:,channel), active, ...
            sqrt(mean(contribution.^2,1)).', sqrt(mean(other.^2,1)).', ...
            'VariableNames', {'Neuron','Weight','Active','ChannelRMS','OtherInputsRMS'});
        connections.RelativeRMS = connections.ChannelRMS ./ max(connections.OtherInputsRMS, realmin);
        clear contribution other;
        baselineRawChannel = rawChannel;
        referenceNormInput = norm(normalized, 'fro');
        referenceFeatures = normalized;
    end
    featureStats = [featureStats; table(degrees(index), mean(rawChannel), std(rawChannel), ...
        min(rawChannel), max(rawChannel), mean(normalizedChannel), std(normalizedChannel), ...
        min(normalizedChannel), max(normalizedChannel), ...
        norm(normalizedChannel-referenceChannel)/sqrt(numel(u)), ...
        norm(normalized-referenceFeatures,'fro')/referenceNormInput, ...
        'VariableNames', {'Degrees','RawMean','RawStd','RawMin','RawMax','NormalizedMean', ...
        'NormalizedStd','NormalizedMin','NormalizedMax','ChannelChangeRMS','AllFeaturesRelativeChange'})]; %#ok<AGROW>
    % Intervention on one channel only, with frozen weights/statistics: causal sensitivity, not retraining.
    normalized(channel,:) = -stats.muX(channel)/stats.sigmaX(channel);
    zeroChannel = predictPNNN(model.network, normalized.', stats, rotation(:));
    if index == 1, zeroReference = zeroChannel; end
    signals = {predicted(1).yvalmod, predicted(2).yvalmod, predicted(3).yvalmod, zeroChannel};
    references = {baseline(1).yvalmod, baseline(2).yvalmod, baseline(3).yvalmod, zeroReference};
    labels = ["GMP-DOMP", "PN-IQ-GMP", "PNNN baseline", "PNNN frozen-weight zero-channel intervention"];
    target = phase*y;
    for family = 1:4
        prediction = double(signals{family});
        reference = double(references{family});
        row = table(degrees(index), labels(family), ...
            norm(prediction-phase*reference)/norm(reference), ...
            nmseComplexDb(target(id), prediction(id)), ...
            nmseComplexDb(target(evaluation), prediction(evaluation)), ...
            nmseComplexDb(target, prediction), ...
            'VariableNames', {'Degrees','Model','EquivarianceError','IdentificationNMSEdB','EvaluationNMSEdB','FullSignalNMSEdB'});
        rotationTable = [rotationTable; row]; %#ok<AGROW>
    end
    fprintf('P%d rotation %d degrees complete.\n', bundle.P, degrees(index));
end
report = struct('modelFile', modelFile, 'measurementFile', measurementFile, ...
    'P', bundle.P, 'sampleCount', numel(u), 'channel', channel, ...
    'featureScale', featureScale, 'featureStats', featureStats, ...
    'connections', connections, 'rotationTable', rotationTable, ...
    'interventionRelativeChange', norm(double(zeroReference)-double(baseline(3).yvalmod))/norm(double(baseline(3).yvalmod)), ...
    'baselineMetrics', bundle.metrics, 'storedSigma', stats.sigmaX(channel));
resultsRoot = fullfile(root, 'results', 'pnnn_phase_diagnostic');
if ~isfolder(resultsRoot), mkdir(resultsRoot); end
report.directory = tempname(resultsRoot);
mkdir(report.directory);
for field = ["featureScale", "featureStats", "connections", "rotationTable"]
    writetable(report.(field), fullfile(report.directory, field + ".csv"));
end
save(fullfile(report.directory, 'diagnosis.mat'), 'report', 'baseline', ...
    'baselineRawChannel', 'referenceChannel', '-v7.3');
disp(connections);
disp(rotationTable);
fprintf('Diagnosis saved: %s\n', report.directory);
end
