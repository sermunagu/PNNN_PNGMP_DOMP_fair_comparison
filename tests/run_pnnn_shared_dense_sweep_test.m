% Test independent sparse targets derived from one immutable dense PNNN.
% The minimal fixture uses zero fine-tuning epochs so it exercises pruning,
% signatures, masks, costs, and predictions without a long training run.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'complexity'));
addpath(fullfile(projectRoot, 'toolbox', 'metrics'));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn', 'pruning'));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

rng(714, 'twister');
n = 48;
inputDimension = 84;
hiddenNeurons = 3;
features = randn(n, inputDimension);
targets = randn(n, 2);
rotation = ones(n, 1);
y = targets(:, 1) + 1j*targets(:, 2);
split.identificationIndices = (1:32).';
split.fullSignalIndices = (1:n).';
normalization = computePNNNNormalization( ...
    features(split.identificationIndices, :), ...
    targets(split.identificationIndices, :));
network = dlnetwork(buildFairPNNNLayers(inputDimension, hiddenNeurons));
parameterCount = countPNNNParameters(inputDimension, hiddenNeurons);
denseFit = struct('network', network, 'normalization', normalization, ...
    'hiddenNeurons', hiddenNeurons, 'bestDenseEpoch', 1, ...
    'trainingTimeSeconds', 0, 'parameterCount', parameterCount);

cfg = getFairDOMPComparisonConfig(projectRoot);
cfg.pnnn.sparseBaseHiddenNeurons = hiddenNeurons;
cfg.pruning.fineTuneEpochs = 0;
cfg.pruning.fineTuneLearnRateDropPeriod = 1;
cfg.training.learnRateDropPeriod = 1;
runtimeConfig = cfg;
denseSource = struct('denseFit', denseFit, ...
    'signature', buildNetworkSignature(denseFit), ...
    'fineTuneEpochs', 0, ...
    'fineTuneBudgetSource', ...
        "historical parameter-matched 344 selection", ...
    'runtimeConfig', runtimeConfig);
targetsToFit = [100 150];
signatureBefore = buildNetworkSignature(denseFit);
points = runPNNNSparseSweep(denseSource, targetsToFit, ...
    features, targets, rotation, y, split, cfg);
signatureAfter = buildNetworkSignature(denseFit);

assert(isequaln(signatureBefore, signatureAfter));
assert(numel(points) == numel(targetsToFit));
for index = 1:numel(points)
    point = points(index);
    assert(point.target == targetsToFit(index));
    assert(point.row.ActualRealParameters == targetsToFit(index));
    assert(isnan(point.row.InternalValidationNMSEdB));
    assert(point.row.ActiveBiases == parameterCount.realBiases);
    assert(point.maskIntegrityAfterPruning.ok);
    assert(point.maskIntegrityAfterFineTune.ok);
    assert(isequaln(point.denseSourceSignature, signatureBefore));
    assert(all(isfinite(point.identificationPrediction)));
    assert(all(isfinite(point.fullSignalPrediction)));
end

checkpointDirectory = tempname;
mkdir(checkpointDirectory);
checkpointCleanup = onCleanup(@() rmdir(checkpointDirectory, 's'));
sweepIdentity = struct('digest', "fixture-sweep");
experimentSignature = struct('schemaVersion', 2, ...
    'digest', "fixture-experiment");
for index = 1:numel(points)
    target = targetsToFit(index);
    [~, ~, filename] = updateSweepCheckpoint(checkpointDirectory, ...
        "pnnn", target, sweepIdentity, experimentSignature, points(index));
    assert(endsWith(string(filename), ...
        compose("pnnn_target_%04d.mat", target)));
end
[loadedFirst, reusedFirst, firstFile] = updateSweepCheckpoint( ...
    checkpointDirectory, "pnnn", targetsToFit(1), ...
    sweepIdentity, experimentSignature);
[loadedSecond, reusedSecond] = updateSweepCheckpoint( ...
    checkpointDirectory, "pnnn", targetsToFit(2), ...
    sweepIdentity, experimentSignature);
assert(reusedFirst && reusedSecond);
assert(loadedFirst.target == targetsToFit(1));
assert(loadedSecond.target == targetsToFit(2));
checkpointArtifact = struct('schemaVersion', -1);
save(firstFile, 'checkpointArtifact');
[~, reusedCorrupt] = updateSweepCheckpoint(checkpointDirectory, ...
    "pnnn", targetsToFit(1), sweepIdentity, experimentSignature);
[~, reusedUnaffected] = updateSweepCheckpoint(checkpointDirectory, ...
    "pnnn", targetsToFit(2), sweepIdentity, experimentSignature);
assert(~reusedCorrupt && reusedUnaffected);
incompatibleIdentity = struct('digest', "different-sweep");
[~, reusedIncompatible] = updateSweepCheckpoint(checkpointDirectory, ...
    "pnnn", targetsToFit(2), incompatibleIdentity, experimentSignature);
assert(~reusedIncompatible);
clear checkpointCleanup;

fprintf('PNNN SHARED DENSE SWEEP TEST: PASS\n');
