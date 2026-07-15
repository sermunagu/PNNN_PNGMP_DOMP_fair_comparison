function study = combineComparisonStudies(linearStudy, pnnnStudy, split)
% combineComparisonStudies - Align the linear and neural study outputs.
% Scientific fitting is already complete; this function only builds the
% shared result, parameter, hyperparameter, and complexity tables.

nIdentification = numel(split.identificationIndices);
nFullSignal = numel(split.fullSignalIndices);
linearResults = addLinearResultColumns( ...
    linearStudy.comparisonResults, nIdentification, nFullSignal);
pnnnResults = pnnnStudy.comparisonResults( ...
    :, linearResults.Properties.VariableNames);
comparisonResults = [linearResults; pnnnResults];

linearFLOPs = addSparseExecutionColumns(linearStudy.complexityFLOPs);
pnnnFLOPs = addSparseExecutionColumns(pnnnStudy.complexityFLOPs);
pnnnFLOPs = pnnnFLOPs(:, linearFLOPs.Properties.VariableNames);
complexityFLOPs = [linearFLOPs; pnnnFLOPs];

parameterSummary = comparisonResults(:, { ...
    'Model','NumRealParameters','ParameterMatchedTarget', ...
    'ParameterDifference','ActualActiveParams','ActiveWeights', ...
    'ActiveBiases','WeightSparsityPercent','FinalFitSamples'});
selectedHyperparameters = comparisonResults(:, { ...
    'Model','SelectedLambda','DOMPSupportSize', ...
    'InternalValidationNMSEdB','BestDenseEpoch', ...
    'BestFineTuneEpoch','SelectionMethod'});

mainModels = ["Independent PN-IQ full"; ...
    "Complex GMP DOMP parameter-matched"; "PNNN N12 sparse"];
mainResults = comparisonResults( ...
    ismember(comparisonResults.Model, mainModels), :);

study = struct( ...
    'comparisonResults', comparisonResults, ...
    'complexityFLOPs', complexityFLOPs, ...
    'parameterSummary', parameterSummary, ...
    'selectedHyperparameters', selectedHyperparameters, ...
    'mainResults', mainResults, ...
    'linear', linearStudy, ...
    'pnnn', pnnnStudy);
end

function results = addLinearResultColumns(results, nIdentification, nFull)
n = height(results);
results.NNSeed = NaN(n, 1);
results.BestDenseEpoch = NaN(n, 1);
results.BestFineTuneEpoch = NaN(n, 1);
results.NNHiddenNeurons = NaN(n, 1);
results.TrainingTimeSeconds = NaN(n, 1);
results.TargetActiveParams = results.ParameterMatchedTarget;
results.ActualActiveParams = results.NumRealParameters;
results.ActiveWeights = results.NumRealParameters;
results.ActiveBiases = zeros(n, 1);
results.WeightSparsityPercent = zeros(n, 1);
results.FinalFitSamples = repmat(nIdentification, n, 1);
results.FullSignalSamples = repmat(nFull, n, 1);
results.NormalizationSamples = NaN(n, 1);
end

function flops = addSparseExecutionColumns(flops)
if ismember('DenseExecutionCoreFLOPsPerSample', ...
        flops.Properties.VariableNames)
    return;
end
flops.DenseExecutionCoreFLOPsPerSample = flops.CoreFLOPsPerSample;
flops.IdealSparseCoreFLOPsPerSample = flops.CoreFLOPsPerSample;
flops.DenseMatrixFLOPsPerSample = flops.FLOPsPerSample;
flops.SparseZeroWeightsSkipped = false(height(flops), 1);
flops.IdealSparseRealMultiplicationsPerSample = ...
    flops.RealMultiplicationsPerSample;
flops.IdealSparseRealAdditionsPerSample = ...
    flops.RealAdditionsPerSample;
flops.IdealSparseCostRequiresSparseKernel = false(height(flops), 1);
end
