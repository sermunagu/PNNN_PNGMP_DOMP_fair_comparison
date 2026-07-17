function fit = refitFairPNNNSparse(denseFinalFit, features, targets, phaseRotation, y, identificationRows, fullSignalRows, targetActiveParams, nnSeed, selectedFineTuneEpochs, cfg)
% refitFairPNNNSparse - Prune and refit the final identification PNNN.
% Biases remain protected, zero weights remain frozen, and fixed-epoch
% fine-tuning uses all identification rows without full-signal selection.

validateattributes(targetActiveParams, {'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(selectedFineTuneEpochs, {'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
if ~isstruct(denseFinalFit) || ...
        ~isfield(denseFinalFit, 'network') || ...
        ~isfield(denseFinalFit, 'normalization') || ...
        ~isa(denseFinalFit.network, 'dlnetwork')
    error('refitFairPNNNSparse:InvalidDenseFit', ...
        'A final dense identification fit is required.');
end
if ~all(ismember(identificationRows, fullSignalRows)) || ...
        ~isequal(fullSignalRows(:), (1:size(features, 1)).')
    error('refitFairPNNNSparse:InvalidRows', ...
        'Identification rows must be contained in the full signal.');
end

timer = tic;
pruning_cfg = cfg.pruning;
pruning_cfg.targetMode = "activeTrainableParams";
pruning_cfg.targetActiveTrainableParams = double(targetActiveParams);
pruning_cfg.fineTuneEpochs = double(selectedFineTuneEpochs);
if pruning_cfg.includeBiases || ~pruning_cfg.freezePruned || ...
        ~pruning_cfg.fineTuneEnabled || ...
        string(pruning_cfg.scope) ~= "global" || ...
        string(pruning_cfg.structureMode) ~= "unstructured"
    error('refitFairPNNNSparse:InvalidPruningContract', ...
        'Final pruning must be global, frozen, unstructured, and bias-safe.');
end

[pruning_state, pruning_stats] = createMagnitudePruningMasks( ...
    denseFinalFit.network, pruning_cfg);
pruned_network = applyLearnableMasks( ...
    denseFinalFit.network, pruning_state.masks);
[integrity_after_pruning, pruning_stats] = ...
    checkPruningMaskIntegrity(pruned_network, pruning_state, ...
    pruning_stats, "after_final_pruning");
if ~integrity_after_pruning.ok
    error('refitFairPNNNSparse:InitialMaskViolation', ...
        'A final pruned weight is nonzero after applying its mask.');
end
assertBiasMasksProtected(pruned_network, pruning_state.masks);

normalization = denseFinalFit.normalization;

features_identification = normalizeFeatures(features(identificationRows, :), normalization);
targets_identification = normalizeTargets(targets(identificationRows, :), normalization);

fine_tune_cfg = struct('training', cfg.training, 'pruning', pruning_cfg);
rng(nnSeed + double(pruning_cfg.fineTuneSeedOffset), 'twister');
[network, fine_tune_info, pruning_stats] = fineTunePrunedNetwork( ...
    pruned_network, features_identification, targets_identification, ...
    [], [], fine_tune_cfg, pruning_state, pruning_stats);
[integrity_after_fine_tune, pruning_stats] = ...
    checkPruningMaskIntegrity(network, pruning_state, pruning_stats, ...
    "after_final_fine_tune");
if ~integrity_after_fine_tune.ok
    error('refitFairPNNNSparse:FineTuneMaskViolation', ...
        'A final pruned weight became nonzero during fine-tuning.');
end

counts = summarizeTrainableParameters( ...
    network, pruning_cfg.includeBiases, pruning_state.masks);
if counts.actualActiveTrainableParams ~= targetActiveParams
    error('refitFairPNNNSparse:TargetMismatch', ...
        'Final pruning produced %d active parameters instead of %d.', ...
        counts.actualActiveTrainableParams, targetActiveParams);
end
if counts.activeBiasParams ~= counts.totalBiasParams || ...
        counts.prunedBiasParams ~= 0
    error('refitFairPNNNSparse:BiasProtectionFailure', ...
        'Every final bias must remain active and protected.');
end

features_full_signal = normalizeFeatures(features(fullSignalRows, :), normalization);
identification_prediction = predictPhaseNorm(network,features_identification, normalization, phaseRotation(identificationRows));
full_signal_prediction = predictPhaseNorm(network, features_full_signal, normalization, phaseRotation(fullSignalRows));

fit = denseFinalFit;
fit.network = network;
fit.targetActiveParams = double(targetActiveParams);
fit.actualActiveParams = counts.actualActiveTrainableParams;
fit.activeWeights = counts.activeWeightParams;
fit.activeBiases = counts.activeBiasParams;
fit.weightSparsityPercent = 100*counts.prunedWeightParams / ...
    max(counts.totalWeightParams, 1);
fit.fineTuneUpdates = fine_tune_info.Updates;
fit.bestFineTuneEpoch = double(selectedFineTuneEpochs);
fit.identificationPrediction = identification_prediction;
fit.fullSignalPrediction = full_signal_prediction;
fit.identificationNMSEdB = nmseComplexDb( ...
    y(identificationRows), identification_prediction);
fit.fullSignalNMSEdB = nmseComplexDb( ...
    y(fullSignalRows), full_signal_prediction);
fit.pruningState = pruning_state;
fit.pruningStats = pruning_stats;
fit.finalFineTuneInfo = fine_tune_info;
fit.maskIntegrityAfterPruning = integrity_after_pruning;
fit.maskIntegrityAfterFineTune = integrity_after_fine_tune;
fit.finalFitSamples = numel(identificationRows);
fit.fullSignalSamples = numel(fullSignalRows);
fit.normalizationSamples = numel(identificationRows);
fit.finalContract = struct( ...
    'fitDomain', "identificationIndices", ...
    'normalizationDomain', "identificationIndices", ...
    'fullSignalUsedForTraining', false, ...
    'fixedDenseEpochs', denseFinalFit.bestDenseEpoch, ...
    'fixedFineTuneEpochs', double(selectedFineTuneEpochs), ...
    'biasesProtected', true, ...
    'prunedWeightsFrozen', true);
fit.trainingTimeSeconds = denseFinalFit.trainingTimeSeconds + toc(timer);
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
        error('refitFairPNNNSparse:BiasMaskViolation', ...
            'Bias masks must contain only true values.');
    end
end
end
