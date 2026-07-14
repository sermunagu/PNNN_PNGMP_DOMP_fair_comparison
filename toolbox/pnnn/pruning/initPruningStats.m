function stats = initPruningStats(pruningCfg)
% initPruningStats - Create the metadata container for PNNN pruning.
%
% This function initializes pruning counters, fine-tuning fields, and mask
% integrity fields before optional magnitude pruning is applied.
%
% Inputs:
%   pruningCfg - Validated pruning configuration struct.
%
% Outputs:
%   stats - Struct saved into model/deploy metadata.

stats = struct();
stats.enabled = pruningCfg.enabled;
stats.targetMode = char(pruningCfg.targetMode);
stats.structureMode = char(pruningCfg.structureMode);
stats.structuredRanking = char(pruningCfg.structuredRanking);
stats.structuredTargetPolicy = char(pruningCfg.structuredTargetPolicy);
stats.hybridExactTarget = logical(pruningCfg.hybridExactTarget);
stats.sparsityTarget = pruningCfg.sparsity;
stats.sparsityActual = 0;
stats.scope = char(pruningCfg.scope);
stats.includeBiases = pruningCfg.includeBiases;
stats.freezePruned = pruningCfg.freezePruned;

stats.totalPodableParams = 0;
stats.totalTrainableParams = 0;
stats.totalPrunableParams = 0;
stats.protectedTrainableParams = 0;
stats.totalWeightParams = 0;
stats.totalBiasParams = 0;
stats.targetActiveTrainableParams = NaN;
stats.targetActivePrunableParams = NaN;
stats.targetActiveParamGap = NaN;
stats.prunedPrunableParams = 0;
stats.remainingPrunableParams = 0;
stats.remainingTotalTrainableParams = 0;
stats.actualActiveTrainableParams = 0;
stats.actualPrunedPrunableParams = 0;
stats.actualPrunableSparsity = 0;
stats.activeWeightParams = 0;
stats.activeBiasParams = 0;
stats.prunedWeightParams = 0;
stats.prunedBiasParams = 0;
stats.totalInputFeatures = NaN;
stats.effectiveInputFeatures = NaN;
stats.activeInputFeatures = NaN;
stats.prunedInputFeatures = NaN;
stats.totalFeatureGroups = NaN;
stats.activeFeatureGroups = NaN;
stats.prunedFeatureGroups = NaN;
stats.prunedFeatureGroupNames = strings(0, 1);

stats.numPrunedParams = 0;
stats.numRemainingParams = 0;
stats.parameterNames = strings(0, 1);
stats.parameterTotal = [];
stats.parameterPruned = [];
stats.parameterRemaining = [];

stats.fineTuneEnabled = pruningCfg.fineTuneEnabled;
stats.fineTuneRun = false;
stats.fineTuneEpochs = 0;
stats.fineTuneInitialLearnRate = pruningCfg.fineTuneInitialLearnRate;
stats.fineTuneBestEpoch = NaN;
stats.fineTuneBestValidationLoss = NaN;
stats.fineTuneFinalValidationLoss = NaN;
stats.fineTuneFinalTrainLoss = NaN;
stats.maskViolationCount = 0;
stats.maskViolationMaxAbs = 0;
stats.maskIntegrityOk = true;
stats.maskIntegrityStage = "not_run";
stats.maskIntegrityChecks = struct( ...
    "stage", {}, ...
    "violationCount", {}, ...
    "violationMaxAbs", {}, ...
    "tolerance", {}, ...
    "ok", {});
end
