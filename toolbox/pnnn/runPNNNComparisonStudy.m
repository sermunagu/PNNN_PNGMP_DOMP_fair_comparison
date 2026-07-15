function study = runPNNNComparisonStudy(x, y, split, cfg, targetActiveParams)
% runPNNNComparisonStudy - Select, refit and evaluate the three PNNN models.

%% Build phase-normalized neural features
% Construct one shared feature matrix for selection, refit and evaluation.
[features, targets, rotation] = buildPhaseNormDataset( ...
    x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
features = features.';
targets = targets.';
rotation = rotation(:);
inputDimension = size(features, 2);
h4 = cfg.pnnn.denseControlHiddenNeurons;
n12 = cfg.pnnn.sparseBaseHiddenNeurons;
seed = cfg.pnnn.nnSeeds(1);

%% Configure the fair training budget
% Scale the validated dense and sparse schedules to the shared split.
budget = scalePNNNTrainingBudget(numel(x), ...
    numel(split.internalTrainIndices), cfg.training, cfg.pruning);
runtimeCfg = cfg;
runtimeCfg.training.maxEpochs = budget.denseMaxEpochs;
runtimeCfg.training.learnRateDropPeriod = budget.denseLearnRateDropPeriod;
runtimeCfg.training.validationPatience = budget.denseValidationPatience;
runtimeCfg.pruning.targetActiveTrainableParams = targetActiveParams;
runtimeCfg.pruning.fineTuneEpochs = budget.fineTuneEpochs;
runtimeCfg.pruning.fineTuneLearnRateDropPeriod = ...
    budget.fineTuneLearnRateDropPeriod;

%% Select hyperparameters using internal validation
% Reuse compatible selections or run the existing validation protocol once.
selection = selectHyperparameters(cfg, split, features, targets, rotation, ...
    y, h4, n12, targetActiveParams, seed, runtimeCfg);
fprintf('PNNN hyperparameter selection reused: %s\n', ...
    upper(string(selection.reused)));
fprintf('Selected epochs H4/N12/fine-tune: %d / %d / %d\n', ...
    selection.h4.bestDenseEpoch, selection.n12Dense.bestDenseEpoch, ...
    selection.n12Sparse.bestFineTuneEpoch);

%% Final refits on the identification set
% Fit dense and sparse networks using all identification samples.
identification = split.identificationIndices;
fullSignal = split.fullSignalIndices;
fprintf('\nFinal H4 refit on all %d identification rows...\n', ...
    numel(identification));
h4Fit = refitFairPNNNDense(features, targets, rotation, y, ...
    identification, fullSignal, h4, seed, selection.h4.bestDenseEpoch, ...
    runtimeCfg.training);
fprintf('Final N12 dense refit on all identification rows...\n');
n12DenseFit = refitFairPNNNDense(features, targets, rotation, y, ...
    identification, fullSignal, n12, seed, ...
    selection.n12Dense.bestDenseEpoch, runtimeCfg.training);
fprintf('Final N12 pruning and fixed-epoch identification fine-tuning...\n');
n12SparseFit = refitFairPNNNSparse(n12DenseFit, features, targets, ...
    rotation, y, identification, fullSignal, targetActiveParams, seed, ...
    selection.n12Sparse.bestFineTuneEpoch, runtimeCfg);

%% Count inference operations
% Report dense execution and ideal sparse execution under one convention.
h4FLOPs = countPNNNFLOPs("PNNN H4 dense", inputDimension, h4, ...
    cfg.pnnn.M, cfg.pnnn.orders);
n12DenseFLOPs = countPNNNFLOPs("PNNN N12 dense", inputDimension, n12, ...
    cfg.pnnn.M, cfg.pnnn.orders);
n12SparseFLOPs = countSparsePNNNFLOPs("PNNN N12 sparse", ...
    inputDimension, n12, cfg.pnnn.M, cfg.pnnn.orders, ...
    n12SparseFit.activeWeights, n12SparseFit.activeBiases);

%% Package PNNN study results
% Return model rows, costs, selections and final fits together.
results = [ ...
    resultRow("PNNN H4 dense", h4Fit, selection.h4, ...
        targetActiveParams, h4FLOPs, cfg.pnnn.actType); ...
    resultRow("PNNN N12 dense", n12DenseFit, selection.n12Dense, ...
        targetActiveParams, n12DenseFLOPs, cfg.pnnn.actType); ...
    resultRow("PNNN N12 sparse", n12SparseFit, selection.n12Sparse, ...
        targetActiveParams, n12SparseFLOPs, cfg.pnnn.actType)];

study = struct( ...
    'comparisonResults', results, ...
    'complexityFLOPs', [h4FLOPs; n12DenseFLOPs; n12SparseFLOPs], ...
    'selection', selection, ...
    'budget', budget, ...
    'runtimeConfig', runtimeCfg, ...
    'inputDimension', inputDimension, ...
    'h4DenseFit', h4Fit, ...
    'n12DenseFit', n12DenseFit, ...
    'n12SparseFit', n12SparseFit);
end

function selection = selectHyperparameters(cfg, split, features, targets, ...
    rotation, y, h4, n12, target, seed, runtimeCfg)
for source = reusableSources(cfg).'
    hyperFile = fullfile(source, 'selected_hyperparameters.csv');
    resultFile = fullfile(source, 'comparison_results.csv');
    splitFile = fullfile(source, 'split_indices.mat');
    if ~(isfile(hyperFile) && isfile(resultFile) && isfile(splitFile))
        continue;
    end
    saved = load(splitFile);
    % TODO: Also match measurement identity/configuration before new captures.
    if sameSelectionRows(saved, split)
        hyper = readtable(hyperFile, TextType='string');
        previous = readtable(resultFile, TextType='string');
        selection = struct('reused', true, 'source', string(source));
        selection.h4 = savedSelection(hyper, previous, "PNNN H4 dense");
        selection.n12Dense = savedSelection( ...
            hyper, previous, "PNNN N12 dense");
        selection.n12Sparse = savedSelection( ...
            hyper, previous, "PNNN N12 sparse");
        return;
    end
end

selection = struct('reused', false, 'source', "new internal selection");
h4Fit = fitFairPNNNDenseValidation(features, targets, rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, h4, seed, runtimeCfg.training);
n12Fit = fitFairPNNNDenseValidation(features, targets, rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, n12, seed, runtimeCfg.training);
sparseFit = pruneAndFineTuneFairPNNN(n12Fit, features, targets, rotation, ...
    y, split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, target, seed, runtimeCfg);
selection.h4 = selectedSummary(h4Fit);
selection.n12Dense = selectedSummary(n12Fit);
selection.n12Sparse = selectedSummary(sparseFit);
end

function sources = reusableSources(cfg)
sources = strings(0, 1);
if isfolder(cfg.resultsRoot)
    entries = dir(cfg.resultsRoot);
    entries = entries([entries.isdir] & ~ismember({entries.name}, {'.','..'}));
    [~, order] = sort([entries.datenum], 'descend');
    if ~isempty(order)
        entries = entries(order);
        sources = string(fullfile( ...
            {entries.folder}, {entries.name})).';
    end
end
sources(end+1, 1) = string(cfg.historicalDisjointResultDirectory);
sources = unique(sources, 'stable');
end

function matches = sameSelectionRows(saved, split)
required = {'internal_train_indices','internal_validation_indices', ...
    'identification_indices'};
matches = all(isfield(saved, required)) && ...
    isequal(saved.internal_train_indices(:), ...
    split.internalTrainIndices(:)) && ...
    isequal(saved.internal_validation_indices(:), ...
    split.internalValidationIndices(:)) && ...
    isequal(saved.identification_indices(:), ...
    split.identificationIndices(:));
end

function summary = savedSelection(hyper, previous, model)
hyperRow = hyper.Model == model;
previousRow = previous.Model == model;
if nnz(hyperRow) ~= 1 || nnz(previousRow) ~= 1
    error('runPNNNComparisonStudy:MissingSelectionRow', ...
        'The reusable artifacts must contain exactly one row for %s.', model);
end
summary.bestDenseEpoch = double(hyper.BestDenseEpoch(hyperRow));
summary.bestFineTuneEpoch = double(hyper.BestFineTuneEpoch(hyperRow));
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
summary.internalTrainNMSEdB = ...
    double(previous.InternalTrainNMSEdB(previousRow));
summary.internalValidationNMSEdB = ...
    double(previous.InternalValidationNMSEdB(previousRow));
end

function summary = selectedSummary(fit)
summary = struct( ...
    'bestDenseEpoch', fit.bestDenseEpoch, ...
    'bestFineTuneEpoch', fit.bestFineTuneEpoch, ...
    'internalTrainNMSEdB', fit.trainNMSEdB, ...
    'internalValidationNMSEdB', fit.validationNMSEdB);
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
end

function row = resultRow(modelName, fit, selection, target, flops, activation)
Model = string(modelName);
SelectionMethod = "internal validation; fixed identification refit";
RegularizationMode = "Not applicable";
NumRealParameters = fit.actualActiveParams;
ParameterMatchedTarget = double(target);
ParameterDifference = NumRealParameters - ParameterMatchedTarget;
InternalTrainNMSEdB = selection.internalTrainNMSEdB;
InternalValidationNMSEdB = selection.internalValidationNMSEdB;
IdentificationNMSEdB = fit.identificationNMSEdB;
FullSignalNMSEdB = fit.fullSignalNMSEdB;
SelectedLambda = NaN;
DOMPSupportSize = NaN;
EffectiveFeatureCount = fit.parameterCount.inputDimension;
PhaseNormalization = true;
IQCoupling = "nonlinear two-output PNNN";
RelativePredictionErrorToComplex = NaN;
FLOPsPerSample = flops.FLOPsPerSample;
AdditionalOperationsPerSample = describeAdditionalOperations(flops, activation);
NNSeed = fit.nnSeed;
BestDenseEpoch = fit.bestDenseEpoch;
BestFineTuneEpoch = fit.bestFineTuneEpoch;
NNHiddenNeurons = fit.hiddenNeurons;
TrainingTimeSeconds = fit.trainingTimeSeconds;
TargetActiveParams = double(target);
ActualActiveParams = fit.actualActiveParams;
ActiveWeights = fit.activeWeights;
ActiveBiases = fit.activeBiases;
WeightSparsityPercent = fit.weightSparsityPercent;
FinalFitSamples = fit.finalFitSamples;
FullSignalSamples = fit.fullSignalSamples;
NormalizationSamples = fit.normalizationSamples;
row = table(Model, SelectionMethod, RegularizationMode, NumRealParameters, ...
    ParameterMatchedTarget, ParameterDifference, InternalTrainNMSEdB, ...
    InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, SelectedLambda, DOMPSupportSize, ...
    EffectiveFeatureCount, PhaseNormalization, IQCoupling, ...
    RelativePredictionErrorToComplex, FLOPsPerSample, ...
    AdditionalOperationsPerSample, NNSeed, BestDenseEpoch, ...
    BestFineTuneEpoch, NNHiddenNeurons, TrainingTimeSeconds, ...
    TargetActiveParams, ActualActiveParams, ActiveWeights, ...
    ActiveBiases, WeightSparsityPercent, FinalFitSamples, ...
    FullSignalSamples, NormalizationSamples);
end
