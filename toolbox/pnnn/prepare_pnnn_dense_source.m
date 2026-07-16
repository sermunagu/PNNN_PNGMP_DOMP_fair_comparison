function study = prepare_pnnn_dense_source( ...
    x, y, split, cfg, targetActiveParameters)
% prepare_pnnn_dense_source - Select and refit the shared dense N12 PNNN.
% Internal validation selects the dense and fine-tuning epochs once.
% The final dense source uses every identification sample.

if ~isfield(cfg, 'experimentSignature')
    error('prepare_pnnn_dense_source:MissingExperimentSignature', ...
        'A measurement/configuration signature is required for safe reuse.');
end

%% Build the phase-normalized neural dataset
% These features are shared by dense selection, final refit, and all targets.
[features, neuralTargets, phaseRotation] = buildPhaseNormDataset( ...
    x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
features = features.';
neuralTargets = neuralTargets.';
phaseRotation = phaseRotation(:);
inputDimension = size(features, 2);
hiddenNeurons = cfg.pnnn.sparseBaseHiddenNeurons;
seed = cfg.pnnn.nnSeeds(1);

%% Scale the established training budget to the current split
budget = scalePNNNTrainingBudget(numel(x), ...
    numel(split.internalTrainIndices), cfg.training, cfg.pruning);
runtimeConfig = cfg;
runtimeConfig.training.maxEpochs = budget.denseMaxEpochs;
runtimeConfig.training.learnRateDropPeriod = ...
    budget.denseLearnRateDropPeriod;
runtimeConfig.training.validationPatience = budget.denseValidationPatience;
runtimeConfig.pruning.targetActiveTrainableParams = targetActiveParameters;
runtimeConfig.pruning.fineTuneEpochs = budget.fineTuneEpochs;
runtimeConfig.pruning.fineTuneLearnRateDropPeriod = ...
    budget.fineTuneLearnRateDropPeriod;

%% Select the dense and fine-tuning epochs on internal validation
selection = selectN12Hyperparameters(cfg, split, features, neuralTargets, ...
    phaseRotation, y, hiddenNeurons, targetActiveParameters, seed, ...
    runtimeConfig);
fprintf('PNNN hyperparameter selection reused: %s\n', ...
    upper(string(selection.reused)));
fprintf('Selected epochs N12/fine-tune: %d / %d\n', ...
    selection.n12Dense.bestDenseEpoch, ...
    selection.n12Sparse.bestFineTuneEpoch);

%% Refit one immutable dense N12 source on identification
fprintf('Final N12 dense refit on all identification rows...\n');
denseFit = refitFairPNNNDense(features, neuralTargets, phaseRotation, y, ...
    split.identificationIndices, split.fullSignalIndices, hiddenNeurons, ...
    seed, selection.n12Dense.bestDenseEpoch, runtimeConfig.training);

study = struct('selection', selection, 'budget', budget, ...
    'runtimeConfig', runtimeConfig, 'inputDimension', inputDimension, ...
    'features', features, 'targets', neuralTargets, ...
    'rotation', phaseRotation, 'n12DenseFit', denseFit);
end

function selection = selectN12Hyperparameters(cfg, split, features, ...
    neuralTargets, rotation, y, hiddenNeurons, target, seed, runtimeConfig)
for source = reusableSources(cfg).'
    hyperFile = fullfile(source, 'selected_hyperparameters.csv');
    resultFile = fullfile(source, 'comparison_results.csv');
    splitFile = fullfile(source, 'split_indices.mat');
    configFile = fullfile(source, 'comparison_config.mat');
    if ~(isfile(hyperFile) && isfile(resultFile) && ...
            isfile(splitFile) && isfile(configFile))
        continue;
    end
    savedSplit = load(splitFile);
    savedSignature = loadOptionalExperimentSignature(configFile);
    previous = readtable(resultFile, TextType='string');
    sparseRow = previous.Model == "PNNN N12 sparse";
    if nnz(sparseRow) ~= 1
        error('prepare_pnnn_dense_source:MissingSelectionRow', ...
            'Reusable results need exactly one PNNN N12 sparse row.');
    end
    savedTarget = double(previous.NumRealParameters(sparseRow));
    if isReusablePNNNSelection(savedSplit, savedSignature, split, ...
            cfg.experimentSignature, savedTarget, target)
        hyper = readtable(hyperFile, TextType='string');
        selection = struct('reused', true, 'source', string(source), ...
            'n12Dense', savedSelection(hyper, previous, "PNNN N12 dense"), ...
            'n12Sparse', savedSelection(hyper, previous, "PNNN N12 sparse"));
        return;
    end
end

denseValidationFit = fitFairPNNNDenseValidation( ...
    features, neuralTargets, rotation, y, split.internalTrainIndices, ...
    split.internalValidationIndices, split.identificationIndices, ...
    hiddenNeurons, seed, runtimeConfig.training);
sparseValidationFit = pruneAndFineTuneFairPNNN( ...
    denseValidationFit, features, neuralTargets, rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, target, seed, runtimeConfig);
selection = struct('reused', false, 'source', "new internal selection", ...
    'n12Dense', selectedSummary(denseValidationFit), ...
    'n12Sparse', selectedSummary(sparseValidationFit));
end

function sources = reusableSources(cfg)
sources = strings(0, 1);
if isfolder(cfg.resultsRoot)
    entries = dir(cfg.resultsRoot);
    entries = entries([entries.isdir] & ...
        ~ismember({entries.name}, {'.','..'}));
    [~, order] = sort([entries.datenum], 'descend');
    if ~isempty(order)
        entries = entries(order);
        sources = string(fullfile({entries.folder}, {entries.name})).';
    end
end
sources(end + 1, 1) = string(cfg.historicalDisjointResultDirectory);
sources = unique(sources, 'stable');
end

function summary = savedSelection(hyper, previous, model)
hyperRow = hyper.Model == model;
resultRow = previous.Model == model;
if nnz(hyperRow) ~= 1 || nnz(resultRow) ~= 1
    error('prepare_pnnn_dense_source:MissingSelectionRow', ...
        'Reusable artifacts need exactly one row for %s.', model);
end
summary.bestDenseEpoch = double(hyper.BestDenseEpoch(hyperRow));
summary.bestFineTuneEpoch = double(hyper.BestFineTuneEpoch(hyperRow));
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
summary.internalTrainNMSEdB = ...
    double(previous.InternalTrainNMSEdB(resultRow));
summary.internalValidationNMSEdB = ...
    double(previous.InternalValidationNMSEdB(resultRow));
end

function summary = selectedSummary(fit)
summary = struct('bestDenseEpoch', fit.bestDenseEpoch, ...
    'bestFineTuneEpoch', fit.bestFineTuneEpoch, ...
    'internalTrainNMSEdB', fit.trainNMSEdB, ...
    'internalValidationNMSEdB', fit.validationNMSEdB);
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
end
