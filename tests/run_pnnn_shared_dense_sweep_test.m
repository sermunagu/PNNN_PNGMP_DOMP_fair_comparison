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
        features, targets, rotation, y, split, cfg);
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
    assert(point.denseSourceDigest == signatureBefore);
    assert(sum(cellfun(@nnz, point.mask)) == targetsToFit(index));
    assert(all(isfinite(point.fullSignalPrediction)));
end
assert(buildNetworkSignature(denseFit) == signatureBefore);

n12Count = countPNNNParameters(84, 12);
n12Target340ActiveWeights = 340 - n12Count.realBiases;
n12Target340Sparsity = 100 * ...
    (n12Count.realWeights - n12Target340ActiveWeights) / ...
    n12Count.realWeights;
assert(abs(n12Target340Sparsity - 68.4108527) < 1e-7);

fprintf('PNNN SHARED DENSE SWEEP TEST: PASS\n');
