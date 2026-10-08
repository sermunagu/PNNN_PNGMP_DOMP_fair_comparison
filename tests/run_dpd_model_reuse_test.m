function run_dpd_model_reuse_test
% Small synthetic recovery/inference test; no measurement data or dense training.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'toolbox')));
rng(821, 'twister');
cfg = getFairDOMPComparisonConfig(root);
cfg.sweep.parameterGrid = 340;
cfg.gmp.blockSize = 64;
cfg.training.verbose = false;
cfg.pruning.fineTuneLearnRateDropPeriod = 1;
n = 384;
u = 0.2*complex(randn(n, 1), randn(n, 1)) + (0.03 + 0.02j);
x = u - mean(u);
y = 0.8*x + 0.1*circshift(x, 2).*abs(circshift(x, -1)).^2;
y = y - mean(y);
split.identificationIndices = (1:320).';
split.fullSignalIndices = (1:n).';
% Fixed synthetic supports exercise recovery without a costly DOMP search.
paths = struct('complex', (1:170).', 'pniq', (1:170).');
[linear, models] = run_linear_sweep(x, y, split, cfg, paths);
assert(all([linear.complexTable.ActualRealParameters; ...
    linear.pniqTable.ActualRealParameters] == 340));
assert(2*numel(models{1}.coefficients) == 340);
assert(numel(models{2}.coefficientsI) + numel(models{2}.coefficientsQ) == 340);

[features, targets, rotation] = buildPhaseNormDataset( ...
    x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
features = features.';
targets = targets.';
normalization = struct('muX', mean(features), 'sigmaX', std(features), ...
    'muY', mean(targets), 'sigmaY', std(targets));
normalization.sigmaX(normalization.sigmaX == 0) = 1;
network = dlnetwork([featureInputLayer(84, Name="input"); ...
    fullyConnectedLayer(12, Name="fc1"); sigmoidLayer(Name="sigmoid1"); ...
    fullyConnectedLayer(2, Name="fcOut")]);
denseFit = struct('network', network, 'normalization', normalization);
denseSource = struct('denseFit', denseFit, 'digest', buildNetworkSignature(denseFit), ...
    'fineTuneEpochs', 1, 'runtimeConfig', struct('training', cfg.training, 'pruning', cfg.pruning));
point = fit_sparse_pnnn_target(denseSource, 340, features, targets, rotation(:), y, split, cfg);
models{3} = point.model;
counts = summarizeTrainableParameters(point.model.network, point.model.masks);
assert(counts.activeWeightParams == 326 && counts.activeBiasParams == 14);
assert(counts.totalWeightParams + counts.totalBiasParams == 1046);
assert(sum(cellfun(@nnz, point.model.masks)) == 340);
assert(buildNetworkSignature(denseFit) == denseSource.digest);
assert(buildNetworkSignature(struct('network', point.model.network, ...
    'normalization', normalization)) ~= denseSource.digest);
bundle = struct('P', 340, 'config', cfg, 'models', {models}, ...
    'sampleRateHz', 491520000, 'metrics', [linear.complexTable; linear.pniqTable; point.row]);

folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
modelFile = fullfile(folder, 'models.mat');
save(modelFile, 'bundle', '-v7.3');
outputFile = fullfile(folder, 'execution.mat');
dpd = applyDPDModels(modelFile, u, outputFile);
expected = {linear.predictions.complexFull, linear.predictions.pniqFull, point.fullSignalPrediction};
for k = 1:3
    assert(~isreal(dpd(k).yvalmod) && all(isfinite(dpd(k).yvalmod)));
    assert(numel(dpd(k).yvalmod) == n);
    assert(norm(dpd(k).yvalmod - expected{k}) / max(norm(expected{k}), eps) < 1e-6);
end
saved = load(outputFile);
assert(isequal(saved.dpd, dpd) && ~saved.metadata.rfValidationPerformed);
try
    applyDPDModels(modelFile, u, outputFile);
    error('Expected overwrite protection.');
catch exception
    assert(strcmp(exception.identifier, 'applyDPDModels:WouldOverwrite'));
end
assert(isequaln(load(outputFile), saved));
% Periodic memory and phase restoration commute with circular shift/rotation.
shifted = applyDPDModels(modelFile, 1j*circshift(u, 7));
for k = 1:2
    reference = 1j*circshift(dpd(k).yvalmod, 7);
    relativeError = norm(shifted(k).yvalmod - reference) / max(norm(reference), eps);
    assert(relativeError < 1e-5);
end
% PNNN must reproduce the legacy computation, including near-zero sigmaX.
% That convention can amplify roundoff: do not assume numerical equivariance.
rotatedInput = 1j*circshift(u, 7);
rotatedInput = rotatedInput - mean(rotatedInput);
[f, ~, r] = buildPhaseNormDataset(rotatedInput, zeros(size(rotatedInput)), ...
    cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
channels = predict(point.model.network, (f.' - normalization.muX) ./ normalization.sigmaX);
channels = channels .* normalization.sigmaY + normalization.muY;
reference = conj(r(:)) .* complex(channels(:,1), channels(:,2));
assert(norm(shifted(3).yvalmod - reference) / max(norm(reference), eps) < 1e-6);
fresh = applyDPDModels(modelFile, 0.1*complex(randn(400,1), randn(400,1)));
for k = 1:3
    assert(numel(fresh(k).yvalmod) == 400 && all(isfinite(fresh(k).yvalmod)));
end
% Small actual DOMP fit and explicit support recovery are identical.
cfg.sweep.parameterGrid = 8;
original = run_linear_sweep(x, y, split, cfg);
recovered = run_linear_sweep(x, y, split, cfg, original.paths);
assert(isequaln(original, recovered));
fprintf('DPD model reuse: PASS (P340, save/load, inference, phase/memory, export protection).\n');
end
