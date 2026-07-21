% Test immutable dense sources, independent sparse targets, and the complete
% reduced PNNN training path without relying on existing checkpoints.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'complexity'));
addpath(fullfile(projectRoot, 'toolbox', 'metrics'));
addpath(fullfile(projectRoot, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));
addpath(fullfile(projectRoot, 'toolbox', 'splits'));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

%% Existing contract: every target starts from one immutable dense source
rng(714, 'twister');
n = 48;
inputDimension = 84;
hiddenNeurons = 3;
features = randn(n, inputDimension);
neuralTargets = randn(n, 2);
phaseRotation = ones(n, 1);
y = neuralTargets(:, 1) + 1j*neuralTargets(:, 2);
split.identificationIndices = (1:32).';
split.fullSignalIndices = (1:n).';
identificationFeatures = features(split.identificationIndices, :);
identificationTargets = neuralTargets(split.identificationIndices, :);
normalization = struct( ...
    'muX', mean(identificationFeatures, 1), ...
    'sigmaX', std(identificationFeatures, 0, 1), ...
    'muY', mean(identificationTargets, 1), ...
    'sigmaY', std(identificationTargets, 0, 1));
normalization.sigmaX(normalization.sigmaX == 0) = 1;
normalization.sigmaY(normalization.sigmaY == 0) = 1;
layers = [
    featureInputLayer(inputDimension, Name="input")
    fullyConnectedLayer(hiddenNeurons, Name="fc1")
    sigmoidLayer(Name="sigmoid1")
    fullyConnectedLayer(2, Name="fcOut")
];
network = dlnetwork(layers);
parameterCount = countPNNNParameters(inputDimension, hiddenNeurons);
denseFit = struct('network', network, 'normalization', normalization);

cfg = getFairDOMPComparisonConfig(projectRoot);
cfg.pnnn.sparseBaseHiddenNeurons = hiddenNeurons;
runtimeConfig = struct( ...
    'training', struct( ...
        'miniBatchSize', cfg.training.miniBatchSize, ...
        'learnRateDropFactor', cfg.training.learnRateDropFactor, ...
        'verbose', cfg.training.verbose), ...
    'pruning', struct( ...
        'fineTuneInitialLearnRate', ...
            cfg.pruning.fineTuneInitialLearnRate, ...
        'fineTuneLearnRateDropPeriod', 1, ...
        'fineTuneSeedOffset', cfg.pruning.fineTuneSeedOffset));
denseSource = struct('denseFit', denseFit, ...
    'digest', buildNetworkSignature(denseFit), ...
    'fineTuneEpochs', 0, ...
    'runtimeConfig', runtimeConfig);
targetsToFit = [100 150];
signatureBefore = buildNetworkSignature(denseFit);
for index = 1:numel(targetsToFit)
    point = fit_sparse_pnnn_target(denseSource, targetsToFit(index), ...
        features, neuralTargets, phaseRotation, y, split, cfg);
    assert(height(point.row) == 1);
    assert(point.row.ActualRealParameters == targetsToFit(index));
    assert(isnan(point.row.InternalValidationNMSEdB));
    assert(point.row.ActiveBiases == parameterCount.realBiases);
    expectedSparsity = 100 * ...
        (parameterCount.realWeights - point.row.ActiveWeights) / ...
        parameterCount.realWeights;
    assert(abs(point.row.WeightSparsityPercent - expectedSparsity) < 1e-12);
    assert(point.row.ActiveWeights + point.row.ActiveBiases == ...
        point.row.ActualRealParameters);
    assert(point.row.FineTuneEpochs == 0);
    assert(isfinite(point.row.MaxAbsRealParameter));
    assert(point.row.MaxAbsRealParameter >= 0);
    assert(point.denseSourceDigest == signatureBefore);
    assert(sum(cellfun(@nnz, point.mask)) == targetsToFit(index));
    assert(all(isfinite(point.fullSignalPrediction)));
end
assert(buildNetworkSignature(denseFit) == signatureBefore);

%% Reduced real execution: selection, dense refit, pruning, and prediction
rng(715, 'twister');
sampleCount = 160;
x = complex(randn(sampleCount, 1), randn(sampleCount, 1))/sqrt(2);
y = 0.8*x + 0.12*circshift(x, 1).*abs(circshift(x, 2)).^2;
reducedCfg = getFairDOMPComparisonConfig(projectRoot);
reducedCfg.identificationFraction = 0.75;
reducedCfg.internalTrainFraction = 0.75;
reducedCfg.pnnn.sparseBaseHiddenNeurons = hiddenNeurons;
reducedCfg.training.miniBatchSize = 32;
reducedCfg.training.historicalMaxEpochs = 2;
reducedCfg.training.historicalLearnRateDropPeriod = 1;
reducedCfg.training.historicalValidationPatience = 2;
reducedCfg.training.verbose = false;
reducedCfg.pruning.historicalFineTuneEpochs = 2;
reducedSplit = buildCommonComparisonSplit(x, y, reducedCfg);
referenceTarget = 200;

[sourceA, realFeaturesA, realTargetsA, realRotationA] = ...
    prepare_pnnn_dense_source(x, y, reducedSplit, reducedCfg, ...
    referenceTarget);
[sourceB, realFeaturesB, realTargetsB, realRotationB] = ...
    prepare_pnnn_dense_source(x, y, reducedSplit, reducedCfg, ...
    referenceTarget);
sourceA.runtimeConfig.training.verbose = reducedCfg.training.verbose;
sourceB.runtimeConfig.training.verbose = reducedCfg.training.verbose;
assert(isequal(size(realFeaturesA), [sampleCount 84]));
assert(isequal(size(realTargetsA), [sampleCount 2]));
assert(numel(realRotationA) == sampleCount);
assert(isequal(realFeaturesA, realFeaturesB));
assert(isequal(realTargetsA, realTargetsB));
assert(isequal(realRotationA, realRotationB));
assert(sourceA.digest == sourceB.digest);
assert(sourceA.bestDenseEpoch == sourceB.bestDenseEpoch);
assert(sourceA.fineTuneEpochs == sourceB.fineTuneEpochs);

finalTarget = 100;
immutableDigest = sourceA.digest;
pointA = fit_sparse_pnnn_target(sourceA, finalTarget, ...
    realFeaturesA, realTargetsA, realRotationA, y, reducedSplit, reducedCfg);
pointB = fit_sparse_pnnn_target(sourceA, finalTarget, ...
    realFeaturesA, realTargetsA, realRotationA, y, reducedSplit, reducedCfg);
assert(pointA.row.ActualRealParameters == finalTarget);
assert(pointA.row.ActiveBiases == ...
    countPNNNParameters(84, hiddenNeurons).realBiases);
assert(isequal(pointA.mask, pointB.mask));
assert(isequal(pointA.fullSignalPrediction, pointB.fullSignalPrediction));
assert(all(isfinite(pointA.fullSignalPrediction)));
assert(buildNetworkSignature(sourceA.denseFit) == immutableDigest);

fineTuneOptions = struct( ...
    'Mode', "fixed-epochs", ...
    'TrainingRows', reducedSplit.identificationIndices, ...
    'ValidationRows', [], ...
    'TargetActiveParameters', finalTarget, ...
    'NNSeed', reducedCfg.pnnn.nnSeed, ...
    'FineTuneEpochs', sourceA.fineTuneEpochs, ...
    'Config', sourceA.runtimeConfig);
sparseResult = pruneAndFineTunePNNN(sourceA.denseFit, ...
    realFeaturesA, realTargetsA, fineTuneOptions);
assert(sparseResult.counts.activeWeightParams + ...
    sparseResult.counts.activeBiasParams == finalTarget);
learnables = sparseResult.network.Learnables;
for row = 1:height(learnables)
    value = numericLearnable(learnables.Value{row});
    mask = logical(sparseResult.masks{row});
    assert(all(value(~mask) == 0, 'all'));
    if lower(string(learnables.Parameter(row))) == "bias"
        assert(all(mask, 'all'));
    end
end

n12Count = countPNNNParameters(84, 12);
n12Target340ActiveWeights = 340 - n12Count.realBiases;
n12Target340Sparsity = 100 * ...
    (n12Count.realWeights - n12Target340ActiveWeights) / ...
    n12Count.realWeights;
assert(abs(n12Target340Sparsity - 68.4108527) < 1e-7);

fprintf('PNNN SHARED DENSE SWEEP TEST: PASS\n');

function value = numericLearnable(value)
if isa(value, 'dlarray')
    value = extractdata(value);
end
value = gather(value);
end
