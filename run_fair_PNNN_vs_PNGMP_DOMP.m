% Script: run_fair_PNNN_vs_PNGMP_DOMP
% Run the corrected 4% identification/full-signal DOMP comparison.
% Hyperparameters remain internal to identification; every final fit uses it.

clearvars;
clc;

project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'toolbox', 'metrics'));
addpath(fullfile(project_root, 'toolbox', 'complexity'));
addpath(fullfile(project_root, 'toolbox', 'domp'));
addpath(fullfile(project_root, 'toolbox', 'splits'));
addpath(fullfile(project_root, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(project_root, 'toolbox', 'pnnn'));
addpath(fullfile(project_root, 'toolbox', 'pnnn', 'pruning'));

cfg = getFairDOMPComparisonConfig(project_root);
if cfg.warmStart.enabled || cfg.warmStart.useLatestDeploy
    error('run_fair_PNNN_vs_PNGMP_DOMP:WarmStartForbidden', ...
        'The fair PNNN comparison must initialize final networks from scratch.');
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
result_directory = fullfile(cfg.resultsRoot, timestamp);
if ~isfolder(result_directory)
    mkdir(result_directory);
end
diary_file = fullfile(result_directory, 'run_log.txt');
diary(diary_file);
diary_cleanup = onCleanup(@() diary('off'));
fprintf('\n=== Full-signal PNNN versus PN-GMP DOMP comparison ===\n');
fprintf('Result directory: %s\n', result_directory);

measurement = load(cfg.measurementFile, 'x', 'y');
if ~isfield(measurement, 'x') || ~isfield(measurement, 'y')
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidMeasurement', ...
        'The configured measurement must contain x and y.');
end
[x_raw, y_raw] = selectXYByMapping( ...
    measurement.x, measurement.y, cfg.mappingMode);
x_raw = x_raw(:);
y_raw = y_raw(:);
if isempty(x_raw) || numel(x_raw) ~= numel(y_raw) || ...
        any(~isfinite(x_raw)) || any(~isfinite(y_raw))
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidSignals', ...
        'Modeled-block X and Y must be aligned finite vectors.');
end
if cfg.pnnn.removeDC
    x = x_raw - mean(x_raw);
    y = y_raw - mean(y_raw);
else
    x = x_raw;
    y = y_raw;
end

split = buildCommonComparisonSplit(x, y, cfg);
identification_indices = split.identificationIndices;
full_signal_indices = split.fullSignalIndices;
internal_train_indices = split.internalTrainIndices;
internal_validation_indices = split.internalValidationIndices;
assert(all(ismember(identification_indices, full_signal_indices)));
assert(numel(full_signal_indices) == numel(x));
assert(numel(unique(identification_indices)) == ...
    numel(identification_indices));
fprintf('Mapping: %s (local modeled-block X/Y convention)\n', ...
    cfg.mappingMode);
fprintf(['Internal train=%d | internal validation=%d | ' ...
    'identification=%d | full signal=%d\n'], ...
    numel(internal_train_indices), numel(internal_validation_indices), ...
    numel(identification_indices), numel(full_signal_indices));
fprintf('Identification is contained in full signal: YES\n');

fprintf('\nFitting six linear DOMP models under the corrected protocol...\n');
linear_study = runPNGMPDOMPStudy(x, y, split, cfg);
linear_results = augmentLinearResults( ...
    linear_study.comparisonResults, numel(identification_indices), ...
    numel(full_signal_indices));
target_row = linear_results.Model == "Independent PN-IQ full";
if nnz(target_row) ~= 1
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingDynamicTarget', ...
        'Independent PN-IQ must define one dynamic parameter target.');
end
target_active_params = linear_results.NumRealParameters(target_row);

[pnnn_features, pnnn_targets, phase_rotation] = ...
    buildPhaseNormDataset(x, y, cfg.pnnn.M, cfg.pnnn.orders, ...
    cfg.pnnn.featMode);
pnnn_features = pnnn_features.';
pnnn_targets = pnnn_targets.';
phase_rotation = phase_rotation(:);
input_dimension = size(pnnn_features, 2);
h4 = cfg.pnnn.denseControlHiddenNeurons;
n12 = cfg.pnnn.sparseBaseHiddenNeurons;
h4_parameter_count = countPNNNParameters(input_dimension, h4);
n12_parameter_count = countPNNNParameters(input_dimension, n12);
assert(h4_parameter_count.realParameters == h4*(input_dimension+3)+2);
assert(n12_parameter_count.realParameters == n12*(input_dimension+3)+2);

budget = scalePNNNTrainingBudget(numel(x), ...
    numel(internal_train_indices), cfg.training, cfg.pruning);
runtime_cfg = cfg;
runtime_cfg.training.maxEpochs = budget.denseMaxEpochs;
runtime_cfg.training.learnRateDropPeriod = ...
    budget.denseLearnRateDropPeriod;
runtime_cfg.training.validationPatience = ...
    budget.denseValidationPatience;
runtime_cfg.pruning.targetActiveTrainableParams = target_active_params;
runtime_cfg.pruning.fineTuneEpochs = budget.fineTuneEpochs;
runtime_cfg.pruning.fineTuneLearnRateDropPeriod = ...
    budget.fineTuneLearnRateDropPeriod;
nn_seed = cfg.pnnn.nnSeeds(1);

selection = loadOrSelectPNNNHyperparameters(cfg, split, ...
    pnnn_features, pnnn_targets, phase_rotation, y, h4, n12, ...
    target_active_params, nn_seed, runtime_cfg);
fprintf('PNNN hyperparameter selection reused: %s\n', ...
    upper(string(selection.reused)));
fprintf('Selected epochs H4/N12/fine-tune: %d / %d / %d\n', ...
    selection.h4.bestDenseEpoch, ...
    selection.n12Dense.bestDenseEpoch, ...
    selection.n12Sparse.bestFineTuneEpoch);

fprintf('\nFinal H4 refit on all %d identification rows...\n', ...
    numel(identification_indices));
h4_fit = refitFairPNNNDense(pnnn_features, pnnn_targets, ...
    phase_rotation, y, identification_indices, full_signal_indices, ...
    h4, nn_seed, selection.h4.bestDenseEpoch, runtime_cfg.training);
fprintf('Final N12 dense refit on all identification rows...\n');
n12_dense_fit = refitFairPNNNDense(pnnn_features, pnnn_targets, ...
    phase_rotation, y, identification_indices, full_signal_indices, ...
    n12, nn_seed, selection.n12Dense.bestDenseEpoch, ...
    runtime_cfg.training);
fprintf('Final N12 pruning and fixed-epoch identification fine-tuning...\n');
n12_sparse_fit = refitFairPNNNSparse(n12_dense_fit, pnnn_features, ...
    pnnn_targets, phase_rotation, y, identification_indices, ...
    full_signal_indices, target_active_params, nn_seed, ...
    selection.n12Sparse.bestFineTuneEpoch, runtime_cfg);

h4_flops = countPNNNFLOPs("PNNN H4 dense", input_dimension, h4, ...
    cfg.pnnn.M, cfg.pnnn.orders);
n12_dense_flops = countPNNNFLOPs("PNNN N12 dense", ...
    input_dimension, n12, cfg.pnnn.M, cfg.pnnn.orders);
n12_sparse_flops = countSparsePNNNFLOPs("PNNN N12 sparse", ...
    input_dimension, n12, cfg.pnnn.M, cfg.pnnn.orders, ...
    n12_sparse_fit.activeWeights, n12_sparse_fit.activeBiases);

pnnn_results = [ ...
    buildPNNNResultRow("PNNN H4 dense", h4_fit, selection.h4, ...
        target_active_params, h4_flops); ...
    buildPNNNResultRow("PNNN N12 dense", n12_dense_fit, ...
        selection.n12Dense, target_active_params, n12_dense_flops); ...
    buildPNNNResultRow("PNNN N12 sparse", n12_sparse_fit, ...
        selection.n12Sparse, target_active_params, n12_sparse_flops)];
pnnn_results = pnnn_results(:, linear_results.Properties.VariableNames);
comparison_results = [linear_results; pnnn_results];

linear_flops = addCommonFLOPColumns(linear_study.complexityFLOPs);
h4_flops = addCommonFLOPColumns(h4_flops);
n12_dense_flops = addCommonFLOPColumns(n12_dense_flops);
n12_sparse_flops = addCommonFLOPColumns(n12_sparse_flops);
h4_flops = h4_flops(:, linear_flops.Properties.VariableNames);
n12_dense_flops = n12_dense_flops(:, linear_flops.Properties.VariableNames);
n12_sparse_flops = n12_sparse_flops(:, linear_flops.Properties.VariableNames);
complexity_flops = [linear_flops; h4_flops; ...
    n12_dense_flops; n12_sparse_flops];

parameter_summary = comparison_results(:, { ...
    'Model','NumRealParameters','ParameterMatchedTarget', ...
    'ParameterDifference','ActualActiveParams','ActiveWeights', ...
    'ActiveBiases','WeightSparsityPercent','FinalFitSamples'});
selected_hyperparameters = comparison_results(:, { ...
    'Model','SelectedLambda','DOMPSupportSize', ...
    'InternalValidationNMSEdB','BestDenseEpoch', ...
    'BestFineTuneEpoch','SelectionMethod'});
writetable(comparison_results, ...
    fullfile(result_directory, 'comparison_results.csv'));
writetable(parameter_summary, ...
    fullfile(result_directory, 'parameter_summary.csv'));
writetable(selected_hyperparameters, ...
    fullfile(result_directory, 'selected_hyperparameters.csv'));
writetable(complexity_flops, ...
    fullfile(result_directory, 'complexity_flops.csv'));
writetable(linear_study.flopConvention, ...
    fullfile(result_directory, 'flop_convention.csv'));
writetable(linear_study.dompHistorySummary, ...
    fullfile(result_directory, 'domp_history_summary.csv'));
writetable(linear_study.regressorStructureSummary, ...
    fullfile(result_directory, 'regressor_structure_summary.csv'));

save(fullfile(result_directory, 'split_indices.mat'), ...
    'identification_indices','full_signal_indices', ...
    'internal_train_indices','internal_validation_indices','split');
supports = linear_study.supports;
domp_histories = linear_study.dompHistories;
save(fullfile(result_directory, 'domp_supports.mat'), ...
    'supports','domp_histories');

identification_predictions = linear_study.identificationPredictions;
identification_predictions.pnnnH4Dense = h4_fit.identificationPrediction;
identification_predictions.pnnnN12Dense = ...
    n12_dense_fit.identificationPrediction;
identification_predictions.pnnnN12Sparse = ...
    n12_sparse_fit.identificationPrediction;
full_signal_predictions = linear_study.fullSignalPredictions;
full_signal_predictions.pnnnH4Dense = h4_fit.fullSignalPrediction;
full_signal_predictions.pnnnN12Dense = ...
    n12_dense_fit.fullSignalPrediction;
full_signal_predictions.pnnnN12Sparse = ...
    n12_sparse_fit.fullSignalPrediction;
target_identification = y(identification_indices);
target_full_signal = y(full_signal_indices);
prediction_metadata = struct( ...
    'identificationSamples', numel(identification_indices), ...
    'fullSignalSamples', numel(full_signal_indices), ...
    'identificationIncludedInFullSignal', true, ...
    'fullSignalIsIndependentHoldout', false);
save(fullfile(result_directory, 'full_signal_predictions.mat'), ...
    'target_identification','target_full_signal', ...
    'identification_indices','full_signal_indices', ...
    'identification_predictions','full_signal_predictions', ...
    'prediction_metadata','-v7.3');

comparison_metadata = struct();
comparison_metadata.measurementName = cfg.measurementName;
comparison_metadata.mappingMode = cfg.mappingMode;
comparison_metadata.sharedDCRemoval = cfg.pnnn.removeDC;
comparison_metadata.inputDimension = input_dimension;
comparison_metadata.identificationSamples = numel(identification_indices);
comparison_metadata.fullSignalSamples = numel(full_signal_indices);
comparison_metadata.identificationIncludedInFullSignal = true;
comparison_metadata.fullSignalIsIndependentHoldout = false;
comparison_metadata.targetActiveParams = target_active_params;
comparison_metadata.pnnnSelection = selection;
comparison_metadata.trainingBudget = budget;
comparison_metadata.linearStudySeconds = linear_study.totalTimeSeconds;
comparison_metadata.equivalenceRelativeError = ...
    linear_study.equivalenceRelativeError;
comparison_metadata.h4Dense = stripPNNNFit(h4_fit);
comparison_metadata.n12Dense = stripPNNNFit(n12_dense_fit);
comparison_metadata.n12Sparse = stripPNNNFit(n12_sparse_fit);
save(fullfile(result_directory, 'comparison_results.mat'), ...
    'comparison_results','comparison_metadata');
save(fullfile(result_directory, 'comparison_config.mat'), ...
    'runtime_cfg','cfg');

createComparisonPlots(comparison_results, result_directory);
main_names = ["Independent PN-IQ full"; ...
    "Complex GMP DOMP parameter-matched"; "PNNN N12 sparse"];
main_rows = ismember(comparison_results.Model, main_names);
main_results = comparison_results(main_rows, :);
sparse_cost = complexity_flops( ...
    complexity_flops.Model == "PNNN N12 sparse", :);
report_file = fullfile(result_directory, 'fair_comparison_report.tex');
writeFairReport(report_file, main_results, split, ...
    linear_study.equivalenceRelativeError, ...
    sparse_cost.DenseMatrixFLOPsPerSample);
if cfg.report.compilePDF
    compileFairReport(report_file, result_directory);
end

disp(' ');
disp('=== Corrected full-signal main comparison ===');
disp(main_results(:, {'Model','NumRealParameters', ...
    'IdentificationNMSEdB','FullSignalNMSEdB','FLOPsPerSample', ...
    'AdditionalOperationsPerSample'}));
disp(' ');
disp('=== All corrected-protocol models ===');
disp(comparison_results(:, {'Model','NumRealParameters', ...
    'IdentificationNMSEdB','FullSignalNMSEdB','FLOPsPerSample'}));
fprintf('Complex/coupled full-signal relative error: %.6e\n', ...
    linear_study.equivalenceRelativeError);
fprintf('All final model fits used %d identification rows.\n', ...
    numel(identification_indices));
fprintf(['Full-signal evaluation used %d rows and includes the ' ...
    'identification rows.\n'], numel(full_signal_indices));
fprintf('Results: %s\n', result_directory);
diary off;

function results = augmentLinearResults(results, nIdentification, nFull)
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

function row = buildPNNNResultRow(modelName, fit, selection, target, flops)
Model = string(modelName);
SelectionMethod = "internal validation; fixed identification refit";
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
AdditionalOperationsPerSample = describeAdditionalOperations(flops);
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
row = table(Model, SelectionMethod, NumRealParameters, ...
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

function selection = loadOrSelectPNNNHyperparameters(cfg, split, ...
    features, targets, rotation, y, h4, n12, target, seed, runtimeCfg)
sources = reusableSelectionSources(cfg);
for source_index = 1:numel(sources)
    source = sources(source_index);
    hyper_file = fullfile(source, 'selected_hyperparameters.csv');
    result_file = fullfile(source, 'comparison_results.csv');
    split_file = fullfile(source, 'split_indices.mat');
    if isfile(hyper_file) && isfile(result_file) && isfile(split_file)
        saved = load(split_file);
    else
        continue;
    end
    if matchesSelectionRows(saved, split)
        hyper = readtable(hyper_file, TextType='string');
        previous = readtable(result_file, TextType='string');
        selection = struct('reused', true, 'source', string(source));
        selection.h4 = selectionRow(hyper, previous, "PNNN H4 dense");
        selection.n12Dense = selectionRow( ...
            hyper, previous, "PNNN N12 dense");
        selection.n12Sparse = selectionRow( ...
            hyper, previous, "PNNN N12 sparse");
        return;
    end
end

selection = struct('reused', false, 'source', "new internal selection");
h4_selected = fitFairPNNNDenseValidation(features, targets, rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, h4, seed, runtimeCfg.training);
n12_selected = fitFairPNNNDenseValidation(features, targets, rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, n12, seed, runtimeCfg.training);
n12_sparse_selected = pruneAndFineTuneFairPNNN(n12_selected, ...
    features, targets, rotation, y, split.internalTrainIndices, ...
    split.internalValidationIndices, split.identificationIndices, ...
    target, seed, runtimeCfg);
selection.h4 = selectedFitSummary(h4_selected);
selection.n12Dense = selectedFitSummary(n12_selected);
selection.n12Sparse = selectedFitSummary(n12_sparse_selected);
end

function sources = reusableSelectionSources(cfg)
directory_sources = strings(0, 1);
if isfolder(cfg.resultsRoot)
    entries = dir(cfg.resultsRoot);
    entries = entries([entries.isdir] & ...
        ~ismember({entries.name}, {'.','..'}));
    if ~isempty(entries)
        [~, order] = sort([entries.datenum], 'descend');
        entries = entries(order);
        directory_sources = strings(numel(entries), 1);
        for index = 1:numel(entries)
            directory_sources(index) = string(fullfile( ...
                entries(index).folder, entries(index).name));
        end
    end
end
sources = [directory_sources; ...
    string(cfg.historicalDisjointResultDirectory)];
sources = unique(sources, 'stable');
end

function matches = matchesSelectionRows(saved, split)
if all(isfield(saved, {'internal_train_indices', ...
        'internal_validation_indices','identification_indices'}))
    matches = isequal(saved.internal_train_indices(:), ...
        split.internalTrainIndices(:)) && ...
        isequal(saved.internal_validation_indices(:), ...
        split.internalValidationIndices(:)) && ...
        isequal(saved.identification_indices(:), ...
        split.identificationIndices(:));
elseif all(isfield(saved, {'TRAIN','VAL','FIT_POOL'}))
    matches = isequal(saved.TRAIN(:), split.internalTrainIndices(:)) && ...
        isequal(saved.VAL(:), split.internalValidationIndices(:)) && ...
        isequal(saved.FIT_POOL(:), split.identificationIndices(:));
else
    matches = false;
end
end

function summary = selectionRow(hyper, previous, model)
hyper_row = hyper.Model == model;
previous_row = previous.Model == model;
if nnz(hyper_row) ~= 1 || nnz(previous_row) ~= 1
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingSelectionRow', ...
        'The reusable selection artifact lacks %s.', model);
end
summary = struct();
summary.bestDenseEpoch = double(hyper.BestDenseEpoch(hyper_row));
summary.bestFineTuneEpoch = double(hyper.BestFineTuneEpoch(hyper_row));
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
if ismember('InternalTrainNMSEdB', previous.Properties.VariableNames)
    summary.internalTrainNMSEdB = ...
        double(previous.InternalTrainNMSEdB(previous_row));
    summary.internalValidationNMSEdB = ...
        double(previous.InternalValidationNMSEdB(previous_row));
else
    summary.internalTrainNMSEdB = ...
        double(previous.TrainNMSEdB(previous_row));
    summary.internalValidationNMSEdB = ...
        double(previous.ValidationNMSEdB(previous_row));
end
end

function summary = selectedFitSummary(fit)
summary = struct('bestDenseEpoch', fit.bestDenseEpoch, ...
    'bestFineTuneEpoch', fit.bestFineTuneEpoch, ...
    'internalTrainNMSEdB', fit.trainNMSEdB, ...
    'internalValidationNMSEdB', fit.validationNMSEdB);
if isnan(summary.bestFineTuneEpoch)
    summary.bestFineTuneEpoch = 0;
end
end

function flops = addCommonFLOPColumns(flops)
if ~ismember('DenseExecutionCoreFLOPsPerSample', ...
        flops.Properties.VariableNames)
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
end

function description = describeAdditionalOperations(row)
parts = strings(0, 1);
if row.NumAbsPerSample > 0
    if row.NumSqrtPerSample == row.NumAbsPerSample
        parts(end+1) = sprintf('%d magnitude (including %d sqrt)', ...
            row.NumAbsPerSample, row.NumSqrtPerSample);
    else
        parts(end+1) = sprintf('%d magnitude', row.NumAbsPerSample);
    end
end
if row.NumSqrtPerSample > 0 && ...
        row.NumSqrtPerSample ~= row.NumAbsPerSample
    parts(end+1) = sprintf('%d sqrt', row.NumSqrtPerSample);
end
if row.NumRealDivisionsPerSample > 0
    parts(end+1) = sprintf('%d division', ...
        row.NumRealDivisionsPerSample);
end
if row.NumELUPerSample > 0
    parts(end+1) = sprintf('%d ELU', row.NumELUPerSample);
end
if row.NumExpWorstCasePerSample > 0
    parts(end+1) = sprintf('up to %d exp', ...
        row.NumExpWorstCasePerSample);
end
if row.PhaseNormalizationIncluded
    parts(end+1) = ...
        "phase normalization/restoration arithmetic in FLOPs";
end
if isempty(parts)
    description = "none";
else
    description = strjoin(parts, ', ');
end
end

function summary = stripPNNNFit(fit)
fields = {'network','normalization','identificationPrediction', ...
    'fullSignalPrediction','pruningState','finalTrainingInfo', ...
    'finalFineTuneInfo'};
available = intersect(fields, fieldnames(fit), 'stable');
summary = rmfield(fit, available);
end

function createComparisonPlots(results, resultDirectory)
labels = string(results.Model);
figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1150 680]);
scatter(results.NumRealParameters, results.FullSignalNMSEdB, 65, 'filled');
grid on;
xlabel('Active real parameters');
ylabel('Full-signal NMSE (dB)');
title('Full-signal NMSE versus active real parameters');
text(results.NumRealParameters, results.FullSignalNMSEdB, ...
    "  " + labels, 'Interpreter', 'none', 'FontSize', 8);
exportgraphics(figure_handle, fullfile(resultDirectory, ...
    'comparison_nmse_parameters.png'), 'Resolution', 160);
close(figure_handle);

marker_sizes = 45 + 100*sqrt(results.NumRealParameters ./ ...
    max(results.NumRealParameters));
figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1150 680]);
scatter(results.FLOPsPerSample, results.FullSignalNMSEdB, ...
    marker_sizes, 'filled');
grid on;
xlabel('FLOPs/sample');
ylabel('Full-signal NMSE (dB)');
title('Full-signal NMSE versus inference FLOPs/sample');
annotations = labels + " (p=" + string(results.NumRealParameters) + ")";
text(results.FLOPsPerSample, results.FullSignalNMSEdB, ...
    "  " + annotations, 'Interpreter', 'none', 'FontSize', 8);
exportgraphics(figure_handle, fullfile(resultDirectory, ...
    'comparison_flops.png'), 'Resolution', 160);
close(figure_handle);
end

function writeFairReport(filename, mainResults, split, ...
    equivalenceError, sparseDenseMatrixFLOPs)
file_id = fopen(filename, 'w');
if file_id < 0
    error('run_fair_PNNN_vs_PNGMP_DOMP:ReportOpenFailed', ...
        'Could not create the LaTeX report.');
end
cleanup = onCleanup(@() fclose(file_id));
fprintf(file_id, '\\documentclass[10pt]{article}\n');
fprintf(file_id, '\\usepackage[margin=1.35cm]{geometry}\n');
fprintf(file_id, '\\usepackage{booktabs,amsmath,array,tabularx}\n');
fprintf(file_id, '\\setlength{\\emergencystretch}{2em}\n');
fprintf(file_id, '\\title{Full-Signal PNNN versus PN-GMP DOMP Comparison}\n');
fprintf(file_id, '\\author{}\\date{}\\begin{document}\\maketitle\n');
fprintf(file_id, ['\\paragraph{Objective and methodology.} Complex GMP, ' ...
    'its exactly equivalent coupled PN-IQ representation, independent ' ...
    'PN-IQ controls, and parameter-matched PNNNs use the same modeled-block ' ...
    'X/Y samples. DOMP is used only during linear identification.\n']);
fprintf(file_id, ['\\paragraph{Fair comparison protocol.} The classic ' ...
    'selector provides %d identification rows (approximately 4\\%% of ' ...
    'the capture). An internal %d/%d split selects regularization and ' ...
    'neural epochs. Every final model is fitted on all identification ' ...
    'rows and evaluated on all %d signal rows. \\textbf{The full-signal ' ...
    'evaluation includes the identification samples and therefore must ' ...
    'not be interpreted as an independent generalization test.}\n'], ...
    numel(split.identificationIndices), ...
    numel(split.internalTrainIndices), ...
    numel(split.internalValidationIndices), ...
    numel(split.fullSignalIndices));
fprintf(file_id, ['\\paragraph{Main results.} NMSE uses ' ...
    '$10\\log_{10}(\\sum|y-\\hat y|^2/\\sum|y|^2)$.\n']);
fprintf(file_id, '\\begin{center}\\small\\begin{tabularx}{\\linewidth}{XrrrX}\n');
fprintf(file_id, '\\toprule\n');
fprintf(file_id, ['Model & Real parameters & Full-signal NMSE (dB) & ' ...
    'FLOPs/sample & Additional operations/sample \\\\ \\midrule\n']);
for index = 1:height(mainResults)
    fprintf(file_id, '%s & %d & %.3f & %d & %s \\\\ \n', ...
        latexEscape(mainResults.Model(index)), ...
        mainResults.NumRealParameters(index), ...
        mainResults.FullSignalNMSEdB(index), ...
        mainResults.FLOPsPerSample(index), ...
        latexEscape(mainResults.AdditionalOperationsPerSample(index)));
end
fprintf(file_id, '\\bottomrule\\end{tabularx}\\end{center}\n');
fprintf(file_id, ['\\paragraph{Computational cost.} FLOPs/sample counts ' ...
    'real additions and multiplications required for one complex output, ' ...
    'including feature generation, coefficients or weights, biases, phase ' ...
    'normalization, and phase restoration. Magnitudes, roots, divisions, ' ...
    'exponentials, and ELUs remain separate. DOMP is not an inference ' ...
    'cost. The sparse-PNNN FLOP count assumes that zero weights are skipped. ' ...
    'An implementation that executes the original full matrices will ' ...
    'require %d FLOPs/sample instead.\n'], sparseDenseMatrixFLOPs);
fprintf(file_id, ['\\paragraph{Interpretation and limitations.} Independent ' ...
    'PN-IQ provides the strongest NMSE in the main comparison; all three ' ...
    'main models use exactly the same number of active real parameters; and ' ...
    'the sparse PNNN uses the fewest counted FLOPs when zero weights are ' ...
    'actually skipped, while adding ELU evaluations. Complex GMP DOMP-100 ' ...
    'and coupled PN-IQ remain numerically equivalent with relative error ' ...
    '$%.3e$. These analytical counts do not establish latency, energy, ' ...
    'or FPGA resource superiority.\n'], equivalenceError);
fprintf(file_id, '\\end{document}\n');
end

function value = latexEscape(value)
value = char(string(value));
value = strrep(value, '\\', '\\textbackslash{}');
value = strrep(value, '_', '\\_');
value = strrep(value, '%', '\\%');
value = strrep(value, '&', '\\&');
end

function compileFairReport(reportFile, outputDirectory)
command = sprintf(['pdflatex -interaction=nonstopmode -halt-on-error ' ...
    '-output-directory="%s" "%s"'], outputDirectory, reportFile);
[status, output] = system(command);
if status == 0
    [status, second_output] = system(command);
    output = [output newline second_output];
end
if status ~= 0
    error('run_fair_PNNN_vs_PNGMP_DOMP:LaTeXCompilationFailed', ...
        'pdflatex failed:%s%s', newline, output);
end
[report_directory, report_stem] = fileparts(reportFile);
temporary_files = {fullfile(report_directory, [report_stem '.aux']), ...
    fullfile(report_directory, [report_stem '.log'])};
for index = 1:numel(temporary_files)
    if isfile(temporary_files{index})
        delete(temporary_files{index});
    end
end
end
