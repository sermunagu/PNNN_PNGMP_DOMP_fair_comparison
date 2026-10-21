function audit = auditPNNNPhaseSource(sweepDirectory, diagnosisDirectory)
% Inspect saved masks and the shared dense network; never train or change them.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
data = load(fullfile(sweepDirectory, 'experiment.mat'));
saved = load(fullfile(sweepDirectory, 'sweep_dense_source.mat'));
source = saved.denseSource;
stats = source.denseFit.normalization;
cfg = data.cfg;
maskTable = table();
for P = cfg.sweep.parameterGrid
    saved = load(fullfile(sweepDirectory, sprintf('models_P%04d.mat', P)), 'bundle');
    m = saved.bundle.models{3};
    rows = m.network.Learnables;
    row = find(rows.Layer == "fc1" & rows.Parameter == "Weights");
    w = double(extractdata(rows.Value{row}));
    maskTable = [maskTable; table(P, nnz(m.masks{row}(:,15)), max(abs(w(:,15))), ...
        'VariableNames', {'P','ActiveChannel15Connections','MaxAbsChannel15Weight'})]; %#ok<AGROW>
end
capture = load(cfg.measurementFile, 'x');
u = double(capture.x(:));
id = data.split.identificationIndices;
ev = setdiff((1:numel(u)).', id);
rotations = table();
for degrees = [0 30 60 90 135 180 270]
    phase = exp(1j*degrees*pi/180);
    x = phase*u;
    if cfg.pnnn.removeDC, x = x - mean(x); end
    [features, ~, rotation] = buildPhaseNormDataset(x, zeros(size(x)), ...
        cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
    normalized = (features.' - stats.muX) ./ stats.sigmaX;
    baseline = predictPNNN(source.denseFit.network, normalized, stats, rotation(:));
    normalized(:,15) = -stats.muX(15)/stats.sigmaX(15);
    ablated = predictPNNN(source.denseFit.network, normalized, stats, rotation(:));
    if degrees == 0, baseReference = baseline; zeroReference = ablated; end
    predictions = {baseline, ablated};
    references = {baseReference, zeroReference};
    labels = ["Dense baseline", "Dense frozen-weight zero-channel intervention"];
    target = phase*data.y;
    for k = 1:2
        prediction = double(predictions{k});
        reference = double(references{k});
        rotations = [rotations; table(degrees, labels(k), ...
            norm(prediction-phase*reference)/norm(reference), ...
            nmseComplexDb(target(id),prediction(id)), nmseComplexDb(target(ev),prediction(ev)), ...
            nmseComplexDb(target,prediction), 'VariableNames', ...
            {'Degrees','Model','EquivarianceError','IdentificationNMSEdB','EvaluationNMSEdB','FullSignalNMSEdB'})]; %#ok<AGROW>
    end
end
learnables = source.denseFit.network.Learnables;
row = learnables.Layer == "fc1" & learnables.Parameter == "Weights";
w = double(extractdata(learnables.Value{row}));
audit = struct('maskTable', maskTable, 'denseRotations', rotations, ...
    'denseChannel15Weights', w(:,15), 'denseEpochs', source.bestDenseEpoch, ...
    'fineTuneEpochs', source.fineTuneEpochs, ...
    'denseInterventionRelativeChange', norm(double(baseReference)-double(zeroReference))/norm(double(baseReference)));
% Auxiliary regressors conj(x) and |x| do not have phase degree one.
saved = load(fullfile(sweepDirectory, 'models_P0310.mat'), 'bundle');
noncanonical = cell(2,1);
for k = 1:2
    m = saved.bundle.models{k};
    degree = arrayfun(@(r) numel(r.X)-numel(r.Xconj), m.regPopulation);
    noncanonical{k} = find(degree ~= 1);
end
audit.nonCovariantRegressorPositions = noncanonical;
writetable(maskTable, fullfile(diagnosisDirectory, 'mask_audit.csv'));
writetable(rotations, fullfile(diagnosisDirectory, 'dense_rotations.csv'));
save(fullfile(diagnosisDirectory, 'source_audit.mat'), 'audit', '-v7.3');
disp(rotations);
fprintf('Dense-channel intervention relative output change: %.12g\n', audit.denseInterventionRelativeChange);
end
