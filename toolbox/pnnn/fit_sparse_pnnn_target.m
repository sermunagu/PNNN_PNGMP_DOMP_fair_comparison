function point = fit_sparse_pnnn_target(denseSource, target, features, neuralTargets, phaseRotation, y, split, cfg)
% fit_sparse_pnnn_target - Fit one independent sparse N12 budget.
% Every budget starts from the same immutable dense source and reuses the
% fine-tuning duration selected once at the reference sparse budget.

target = double(target);
referenceDigest = denseSource.digest;

%% 1. Verify that the shared dense source has not changed
if buildNetworkSignature(denseSource.denseFit) ~= referenceDigest
    error('fit_sparse_pnnn_target:DenseSourceChanged', ...
        'Every sparse target must start from the same immutable dense source.');
end
denseCounts = summarizeTrainableParameters(denseSource.denseFit.network);

%% 2. Prune independently from the original dense source
fineTuneOptions = struct( ...
    'Mode', "fixed-epochs", ...
    'TrainingRows', split.identificationIndices, ...
    'ValidationRows', [], ...
    'TargetActiveParameters', target, ...
    'NNSeed', cfg.pnnn.nnSeed, ...
    'FineTuneEpochs', denseSource.fineTuneEpochs, ...
    'Config', denseSource.runtimeConfig);

t = tic;
sparseFit = pruneAndFineTunePNNN(denseSource.denseFit, features, neuralTargets, fineTuneOptions);
fprintf('Fine-tuning: %.1f s\n', toc(t));

if buildNetworkSignature(denseSource.denseFit) ~= referenceDigest
    error('fit_sparse_pnnn_target:DenseSourceChanged', ...
        'A sparse fit modified the immutable dense source.');
end

%% 3. Predict the identification subset and complete signal
normalization = denseSource.denseFit.normalization;

identificationRows = split.identificationIndices;
fullSignalRows = split.fullSignalIndices;

featuresIdentification = (features(identificationRows, :) - normalization.muX) ./ normalization.sigmaX;
featuresFullSignal = (features(fullSignalRows, :) - normalization.muX) ./ normalization.sigmaX;

identificationPrediction = predictPhaseNorm(sparseFit.network, featuresIdentification, normalization, phaseRotation(identificationRows));
fullSignalPrediction = predictPhaseNorm(sparseFit.network, featuresFullSignal, normalization, phaseRotation(fullSignalRows));

%% 4. Calculate error, complexity, sparsity, and parameter range
counts = sparseFit.counts;

identificationNMSEdB = nmseComplexDb( y(identificationRows), identificationPrediction);
fullSignalNMSEdB = nmseComplexDb(y(fullSignalRows), fullSignalPrediction);

flops = countSparsePNNNFLOPs(cfg.names.pnnn, ...
    size(features, 2), cfg.pnnn.sparseBaseHiddenNeurons, ...
    cfg.pnnn.M, cfg.pnnn.orders, counts.activeWeightParams, ...
    counts.activeBiasParams);

weightSparsityPercent = 100 * ...
    (denseCounts.totalWeightParams - counts.activeWeightParams) / ...
    denseCounts.totalWeightParams;

%% 5. Build the unchanged public result row and checkpoint point
Model = cfg.names.pnnn;
TargetRealParameters = target;
ActualRealParameters = double(counts.activeWeightParams + counts.activeBiasParams);
SelectedLambda = NaN;
InternalValidationNMSEdB = NaN;
IdentificationNMSEdB = double(identificationNMSEdB);
FullSignalNMSEdB = double(fullSignalNMSEdB);
FLOPsPerSample = double(flops.FLOPsPerSample);
ActiveWeights = double(counts.activeWeightParams);
ActiveBiases = double(counts.activeBiasParams);
WeightSparsityPercent = double(weightSparsityPercent);
FineTuneEpochs = double(denseSource.fineTuneEpochs);
MaxAbsRealParameter = double(counts.maxAbsRealParameter);
row = table(Model, TargetRealParameters, ActualRealParameters, ...
    SelectedLambda, InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, ActiveWeights, ActiveBiases, ...
    WeightSparsityPercent, FineTuneEpochs, MaxAbsRealParameter);

point = struct('row', row, 'mask', {sparseFit.masks}, ...
    'fullSignalPrediction', fullSignalPrediction, ...
    'denseSourceDigest', referenceDigest);
end
