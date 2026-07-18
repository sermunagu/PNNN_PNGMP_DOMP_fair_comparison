function bestFineTuneEpoch = pruneAndFineTuneFairPNNN( ...
    denseFit, features, targets, trainingRows, validationRows, ...
    targetActiveParameters, nnSeed, cfg)
% pruneAndFineTuneFairPNNN - Select the shared sparse fine-tuning epoch.

masks = createMagnitudePruningMasks( ...
    denseFit.network, targetActiveParameters);
network = applyLearnableMasks(denseFit.network, masks);
checkPruningMaskIntegrity(network, masks);

normalization = denseFit.normalization;
featuresTraining = normalizeFeatures( ...
    features(trainingRows, :), normalization);
featuresValidation = normalizeFeatures( ...
    features(validationRows, :), normalization);
targetsTraining = normalizeTargets(targets(trainingRows, :), normalization);
targetsValidation = normalizeTargets( ...
    targets(validationRows, :), normalization);

rng(nnSeed + double(cfg.pruning.fineTuneSeedOffset), 'twister');
[network, bestFineTuneEpoch] = fineTunePrunedNetwork( ...
    network, featuresTraining, targetsTraining, ...
    featuresValidation, targetsValidation, cfg, masks);
checkPruningMaskIntegrity(network, masks);

counts = summarizeTrainableParameters(network, masks);
if counts.activeWeightParams + counts.activeBiasParams ~= ...
        targetActiveParameters
    error('pruneAndFineTuneFairPNNN:TargetMismatch', ...
        'Pruning did not produce the requested active parameter count.');
end
if counts.activeBiasParams ~= counts.totalBiasParams
    error('pruneAndFineTuneFairPNNN:BiasProtectionFailure', ...
        'Every bias must remain active and protected.');
end
end

function values = normalizeFeatures(values, stats)
values = (values - stats.muX) ./ stats.sigmaX;
end

function values = normalizeTargets(values, stats)
values = (values - stats.muY) ./ stats.sigmaY;
end
