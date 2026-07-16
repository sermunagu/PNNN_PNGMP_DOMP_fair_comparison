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
    assert(height(point.row) == 1);
    assert(isscalar(point.row.ArtifactFile));
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
nameTargets = [300 344 380];
artifactNames = strings(size(nameTargets));
nameDirectory = fullfile(checkpointDirectory, 'artifact_names');
for index = 1:numel(nameTargets)
    [~, ~, filename] = updateSweepCheckpoint(nameDirectory, ...
        "pnnn", nameTargets(index), sweepIdentity, experimentSignature, ...
        struct('target', nameTargets(index)));
    [~, base, extension] = fileparts(filename);
    artifactNames(index) = string(base) + string(extension);
    assert(isscalar(artifactNames(index)));
end
assert(numel(unique(artifactNames)) == numel(nameTargets));
assert(isequal(artifactNames, compose( ...
    "pnnn_target_%04d.mat", nameTargets)));

reusable = true;
reusedPoint = points(1);
reusedPoint.row.ArtifactFile = ["pnnn_target_0300", ".mat"];
if reusable
    legacyArtifactFile = reusedPoint.row.ArtifactFile;
    assert(isequal(size(legacyArtifactFile), [1 2]));
    assert(legacyArtifactFile(1) + legacyArtifactFile(2) == ...
        artifactNames(1));
    reusedPoint.row.ArtifactFile = artifactNames(1);
end
assert(isscalar(reusedPoint.row.ArtifactFile));
assert(reusedPoint.row.ArtifactFile == "pnnn_target_0300.mat");

schema = cfg.sweep.resultSchema;
pnnnRows = table('Size', [0 numel(schema.names)], ...
    'VariableTypes', schema.types, 'VariableNames', schema.names);
pointRow = points(1).row;
pointRow.ArtifactFile = artifactNames(1);
assert(height(pointRow) == 1 && width(pointRow) == width(pnnnRows));
assert(isequal(pointRow.Properties.VariableNames, ...
    pnnnRows.Properties.VariableNames));
for index = 1:width(pointRow)
    name = schema.names{index};
    assert(string(class(pointRow.(name))) == schema.types(index));
end
pnnnRows(1, :) = pointRow;
assert(isscalar(pnnnRows.ArtifactFile));

checkpointDirectories = { ...
    fullfile(checkpointDirectory, 'char_input'), ...
    string(fullfile(checkpointDirectory, 'string_input'))};
artifactKinds = ["linear", "dense", "pnnn"];
artifactTargets = [NaN NaN targetsToFit(1)];
expectedFiles = sort(["linear_sweep.mat", ...
    "sweep_dense_source.mat", "pnnn_target_0100.mat"]);
for pathIndex = 1:numel(checkpointDirectories)
    for artifactIndex = 1:numel(artifactKinds)
        payload = struct('kind', artifactKinds(artifactIndex), ...
            'pathIndex', pathIndex);
        [~, ~, filename] = updateSweepCheckpoint( ...
            checkpointDirectories{pathIndex}, artifactKinds(artifactIndex), ...
            artifactTargets(artifactIndex), sweepIdentity, ...
            experimentSignature, payload);
        assert(isstring(filename) && isscalar(filename));
        assert(isfile(char(filename)));
        [loaded, reusable] = updateSweepCheckpoint( ...
            checkpointDirectories{pathIndex}, artifactKinds(artifactIndex), ...
            artifactTargets(artifactIndex), sweepIdentity, ...
            experimentSignature);
        assert(reusable && isequaln(loaded, payload));
    end
    savedFiles = dir(fullfile( ...
        char(checkpointDirectories{pathIndex}), '*.mat'));
    assert(isequal(sort(string({savedFiles.name})), expectedFiles));
end
invalidDirectory = [string(fullfile(checkpointDirectory, 'invalid_a')), ...
    string(fullfile(checkpointDirectory, 'invalid_b'))];
invalidPathRejected = false;
try
    updateSweepCheckpoint(invalidDirectory, "linear", [], ...
        sweepIdentity, experimentSignature, struct('value', 1));
catch exception
    invalidPathRejected = strcmp(exception.identifier, ...
        'updateSweepCheckpoint:InvalidPath');
end
assert(invalidPathRejected);

pointDirectory = fullfile(checkpointDirectory, 'pnnn_points');
for index = 1:numel(points)
    target = targetsToFit(index);
    [~, ~, filename] = updateSweepCheckpoint(pointDirectory, ...
        "pnnn", target, sweepIdentity, experimentSignature, points(index));
    assert(endsWith(string(filename), ...
        compose("pnnn_target_%04d.mat", target)));
end
[loadedFirst, reusedFirst, firstFile] = updateSweepCheckpoint( ...
    pointDirectory, "pnnn", targetsToFit(1), ...
    sweepIdentity, experimentSignature);
[loadedSecond, reusedSecond] = updateSweepCheckpoint( ...
    pointDirectory, "pnnn", targetsToFit(2), ...
    sweepIdentity, experimentSignature);
assert(reusedFirst && reusedSecond);
assert(loadedFirst.target == targetsToFit(1));
assert(loadedSecond.target == targetsToFit(2));
checkpointArtifact = struct('schemaVersion', -1);
save(char(firstFile), 'checkpointArtifact');
[~, reusedCorrupt] = updateSweepCheckpoint(pointDirectory, ...
    "pnnn", targetsToFit(1), sweepIdentity, experimentSignature);
[~, reusedUnaffected] = updateSweepCheckpoint(pointDirectory, ...
    "pnnn", targetsToFit(2), sweepIdentity, experimentSignature);
assert(~reusedCorrupt && reusedUnaffected);
incompatibleIdentity = struct('digest', "different-sweep");
[~, reusedIncompatible] = updateSweepCheckpoint(pointDirectory, ...
    "pnnn", targetsToFit(2), incompatibleIdentity, experimentSignature);
assert(~reusedIncompatible);
clear checkpointCleanup;

fprintf('PNNN SHARED DENSE SWEEP TEST: PASS\n');
