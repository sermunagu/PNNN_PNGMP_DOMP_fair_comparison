function point = fit_sparse_pnnn_target(denseSource, target, ...
    features, targets, rotation, y, split, cfg)
% fit_sparse_pnnn_target - Fit independent sparse N12 points.
% Every target starts from the same immutable dense source; biases stay
% protected and the shared fine-tuning budget is reused without reselection.

target = double(target);
referenceDigest = denseSource.digest;
if buildNetworkSignature(denseSource.denseFit) ~= referenceDigest
    error('fit_sparse_pnnn_target:DenseSourceChanged', ...
        'Every sparse target must start from the same immutable dense source.');
end

denseCounts = summarizeTrainableParameters(denseSource.denseFit.network);

%% Prune, refit, predict, and count one independent target
% refitFairPNNNSparse contains the mathematical pruning and frozen-mask update.
fit = refitFairPNNNSparse(denseSource.denseFit, features, targets, ...
    rotation, y, split.identificationIndices, split.fullSignalIndices, ...
    target, cfg.pnnn.nnSeed, denseSource.fineTuneEpochs, ...
    denseSource.runtimeConfig);
if buildNetworkSignature(denseSource.denseFit) ~= referenceDigest
    error('fit_sparse_pnnn_target:DenseSourceChanged', ...
        'A sparse fit modified the immutable dense source.');
end

flops = countSparsePNNNFLOPs("Sparse PNNN N12", ...
    size(features, 2), cfg.pnnn.sparseBaseHiddenNeurons, ...
    cfg.pnnn.M, cfg.pnnn.orders, fit.activeWeights, fit.activeBiases);
weightSparsityPercent = 100 * ...
    (denseCounts.totalWeightParams - fit.activeWeights) / ...
    denseCounts.totalWeightParams;

Model = "Sparse PNNN N12";
TargetRealParameters = target;
ActualRealParameters = double(fit.activeWeights + fit.activeBiases);
SelectedLambda = NaN;
InternalValidationNMSEdB = NaN;
IdentificationNMSEdB = double(fit.identificationNMSEdB);
FullSignalNMSEdB = double(fit.fullSignalNMSEdB);
FLOPsPerSample = double(flops.FLOPsPerSample);
ActiveWeights = double(fit.activeWeights);
ActiveBiases = double(fit.activeBiases);
WeightSparsityPercent = double(weightSparsityPercent);
FineTuneEpochs = double(denseSource.fineTuneEpochs);
row = table(Model, TargetRealParameters, ActualRealParameters, ...
    SelectedLambda, InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, ActiveWeights, ActiveBiases, ...
    WeightSparsityPercent, FineTuneEpochs);

point = struct('row', row, 'mask', {fit.masks}, ...
    'fullSignalPrediction', fit.fullSignalPrediction, ...
    'denseSourceDigest', referenceDigest);
end
