function fit = pruneAndFineTuneFairPNNN( ...
    denseFit, features, targets, phaseRotation, y, trainingRows, ...
    validationRows, fitPoolRows, targetActiveParams, nnSeed, cfg)
% pruneAndFineTuneFairPNNN - Prune N12 and fine-tune with frozen masks.
% Global magnitude masks protect every bias and target an exact active
% trainable-parameter count derived by the caller from Independent PN-IQ.

validateattributes(targetActiveParams, {'numeric'}, ...
    {'scalar','integer','positive','finite'});
if ~isstruct(denseFit) || ~isfield(denseFit, 'network') || ...
        ~isfield(denseFit, 'normalization') || ...
        ~isa(denseFit.network, 'dlnetwork')
    error('pruneAndFineTuneFairPNNN:InvalidDenseFit', ...
        'denseFit must contain the trained dense dlnetwork and normalization.');
end

pruning_cfg = cfg.pruning;
pruning_cfg.targetMode = "activeTrainableParams";
pruning_cfg.targetActiveTrainableParams = double(targetActiveParams);
if pruning_cfg.includeBiases || ~pruning_cfg.freezePruned || ...
        ~pruning_cfg.fineTuneEnabled || ...
        string(pruning_cfg.scope) ~= "global" || ...
        string(pruning_cfg.structureMode) ~= "unstructured"
    error('pruneAndFineTuneFairPNNN:InvalidPruningContract', ...
        'The fair comparison requires global frozen unstructured weight pruning with protected biases.');
end

timer = tic;
[pruning_state, pruning_stats] = createMagnitudePruningMasks( ...
    denseFit.network, pruning_cfg);
pruned_network = applyLearnableMasks( ...
    denseFit.network, pruning_state.masks);
[integrity_after_pruning, pruning_stats] = checkPruningMaskIntegrity( ...
    pruned_network, pruning_state, pruning_stats, "after_pruning");
if ~integrity_after_pruning.ok
    error('pruneAndFineTuneFairPNNN:InitialMaskViolation', ...
        'A pruned weight is nonzero immediately after applying the mask.');
end
assertBiasMasksProtected(pruned_network, pruning_state.masks);

normalization = denseFit.normalization;
features_training = normalizeFeatures(features(trainingRows, :), ...
    normalization);
features_validation = normalizeFeatures(features(validationRows, :), ...
    normalization);
targets_training = normalizeTargets(targets(trainingRows, :), ...
    normalization);
targets_validation = normalizeTargets(targets(validationRows, :), ...
    normalization);

fine_tune_cfg = struct();
fine_tune_cfg.training = cfg.training;
fine_tune_cfg.pruning = pruning_cfg;
rng(nnSeed + double(pruning_cfg.fineTuneSeedOffset), 'twister');
[network, fine_tune_info, pruning_stats] = fineTunePrunedNetwork( ...
    pruned_network, features_training, targets_training, ...
    features_validation, targets_validation, fine_tune_cfg, ...
    pruning_state, pruning_stats);
[integrity_after_fine_tune, pruning_stats] = checkPruningMaskIntegrity( ...
    network, pruning_state, pruning_stats, "after_fine_tune");
if ~integrity_after_fine_tune.ok
    error('pruneAndFineTuneFairPNNN:FineTuneMaskViolation', ...
        'A pruned weight became nonzero during fine-tuning.');
end

counts = summarizeTrainableParameters( ...
    network, pruning_cfg.includeBiases, pruning_state.masks);
if counts.actualActiveTrainableParams ~= targetActiveParams
    error('pruneAndFineTuneFairPNNN:TargetMismatch', ...
        'Pruning produced %d active parameters instead of %d.', ...
        counts.actualActiveTrainableParams, targetActiveParams);
end
if counts.activeBiasParams ~= counts.totalBiasParams || ...
        counts.prunedBiasParams ~= 0
    error('pruneAndFineTuneFairPNNN:BiasProtectionFailure', ...
        'Every bias must remain active and protected.');
end

fit = denseFit;
fit.network = network;
fit.targetActiveParams = double(targetActiveParams);
fit.actualActiveParams = counts.actualActiveTrainableParams;
fit.activeWeights = counts.activeWeightParams;
fit.activeBiases = counts.activeBiasParams;
fit.weightSparsityPercent = 100*counts.prunedWeightParams / ...
    max(counts.totalWeightParams, 1);
fit.fineTuneUpdates = fine_tune_info.Updates;
fit.bestFineTuneEpoch = fine_tune_info.BestEpoch;
fit.bestEpoch = fine_tune_info.BestEpoch;
fit.trainNMSEdB = predictNMSE(network, features, phaseRotation, y, ...
    trainingRows, normalization);
fit.validationNMSEdB = predictNMSE(network, features, phaseRotation, y, ...
    validationRows, normalization);
fit.fitPoolNMSEdB = predictNMSE(network, features, phaseRotation, y, ...
    fitPoolRows, normalization);
fit.pruningState = pruning_state;
fit.pruningStats = pruning_stats;
fit.fineTuneInfo = fine_tune_info;
fit.maskIntegrityAfterPruning = integrity_after_pruning;
fit.maskIntegrityAfterFineTune = integrity_after_fine_tune;
fit.trainingTimeSeconds = denseFit.trainingTimeSeconds + toc(timer);
fit.selectionContract.testRowsUsed = false;
fit.selectionContract.pruningRankingDomain = "trained dense weights";
fit.selectionContract.fineTuneDomain = "TRAIN";
fit.selectionContract.fineTuneCheckpointDomain = "VAL";
end

function value = predictNMSE(network, features, rotation, y, rows, stats)
features_normalized = normalizeFeatures(features(rows, :), stats);
prediction = predictPhaseNorm(network, features_normalized, stats, ...
    rotation(rows));
value = nmseComplexDb(y(rows), prediction);
end

function values = normalizeFeatures(values, stats)
values = (values - stats.muX) ./ stats.sigmaX;
end

function values = normalizeTargets(values, stats)
values = (values - stats.muY) ./ stats.sigmaY;
end

function assertBiasMasksProtected(network, masks)
learnables = network.Learnables;
for row = 1:height(learnables)
    if lower(string(learnables.Parameter(row))) == "bias" && ...
            ~all(logical(masks{row}), 'all')
        error('pruneAndFineTuneFairPNNN:BiasMaskViolation', ...
            'Bias masks must contain only true values.');
    end
end
end
