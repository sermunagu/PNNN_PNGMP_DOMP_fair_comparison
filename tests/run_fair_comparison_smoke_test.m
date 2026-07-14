% Script: run_fair_comparison_smoke_test
% Verify the 4% identification/full-signal contract, DOMP-only paths,
% fixed final PNNN refits, phase restoration, and shared FLOP scope.

clearvars;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(project_root, 'toolbox', 'domp'));
addpath(fullfile(project_root, 'toolbox', 'complexity'));
addpath(fullfile(project_root, 'toolbox', 'metrics'));
addpath(fullfile(project_root, 'toolbox', 'splits'));
addpath(fullfile(project_root, 'toolbox', 'pnnn'));
addpath(fullfile(project_root, 'toolbox', 'pnnn', 'pruning'));

rng(712, 'twister');
n_samples = 768;
x = (randn(n_samples, 1) + 1j*randn(n_samples, 1))/sqrt(2);
x(1:97:end) = 0;
y = 0.82*x + 0.13*circshift(x, 1).*abs(circshift(x, 2)).^2 + ...
    0.04*conj(circshift(x, 3));

cfg = getFairDOMPComparisonConfig(project_root);
smoke_cfg = cfg.gmp;
smoke_cfg.maxPopulation = 12;
smoke_cfg.blockSize = 128;
smoke_cfg.selectionMethod = 'DOMP';
smoke_cfg.lambda1 = 1e-3;
smoke_cfg.lambda2 = 1e-4;

split = buildCommonComparisonSplit(x, y, cfg);
split_repeat = buildCommonComparisonSplit(x, y, cfg);
assert(isequal(split.internalTrainIndices, ...
    split_repeat.internalTrainIndices));
assert(isequal(split.internalValidationIndices, ...
    split_repeat.internalValidationIndices));
assert(isequal(split.identificationIndices, ...
    split_repeat.identificationIndices));
assert(isequal(split.fullSignalIndices, split_repeat.fullSignalIndices));
identification_indices = split.identificationIndices;
full_signal_indices = split.fullSignalIndices;
assert(abs(numel(identification_indices) - ...
    floor(cfg.identificationFraction*n_samples)) <= 2);
assert(all(ismember(identification_indices, full_signal_indices)));
assert(numel(full_signal_indices) == n_samples);
assert(numel(unique(identification_indices)) == ...
    numel(identification_indices));
assert(isempty(intersect(split.internalTrainIndices, ...
    split.internalValidationIndices)));
assert(isequal(sort([split.internalTrainIndices; ...
    split.internalValidationIndices]), sort(identification_indices)));

% Verify the phase rotation and the exact zero-input convention.
phase_rotation = computePhaseNormGMPRotation(x, full_signal_indices);
x_validation = x(full_signal_indices);
zero_input = abs(x_validation) == 0;
nonzero_input = ~zero_input;
expected_rotation = conj(x_validation(nonzero_input))./ ...
    abs(x_validation(nonzero_input));
assert(iscolumn(phase_rotation) && all(isfinite(phase_rotation)));
assert(all(phase_rotation(zero_input) == 1));
assert(norm(phase_rotation(nonzero_input) - expected_rotation, Inf) <= 1e-13);
rotated_input = phase_rotation.*x_validation;
input_scale = max(1, norm(abs(x_validation), Inf));
assert(norm(real(rotated_input) - abs(x_validation), Inf) ...
    <= 1e-12*input_scale);
assert(norm(imag(rotated_input), Inf) <= 1e-12*input_scale);

% Build one population and obtain the fixture complex support.
gmp_manager = GMP_createRegressorManager(x, y, smoke_cfg);
[~, complex_fit] = GMP_blockFitEvaluate( ...
    x, y, identification_indices, full_signal_indices, ...
    gmp_manager, smoke_cfg, 'Smoke complex GMP');
n_population = numel(gmp_manager.regPopulation);
support_complex = complex_fit.support(:);
assert(n_population == 673);
assert(numel(support_complex) == 12);
assert(strcmpi(complex_fit.selectionMethod, 'DOMP'));

[~, structure_analysis] = analyzeRegressorStructure( ...
    gmp_manager.regPopulation, support_complex);
descriptors = structure_analysis.descriptors;
assert(numel(descriptors) == 673);
assert(nnz([descriptors.canonicalGMP]) == 671);
auxiliary_indices = find(~[descriptors.canonicalGMP]).';
auxiliary_types = string({descriptors(auxiliary_indices).auxiliaryType});
assert(numel(auxiliary_indices) == 2);
assert(isequal(sort(auxiliary_types(:)), sort(["conjugate"; "envelope"])));

all_population = (1:n_population).';
U_identification = buildGMPRegressorRows( ...
    x, identification_indices, gmp_manager, all_population);
phase_rotation_identification = computePhaseNormGMPRotation( ...
    x, identification_indices);
U_phase_normalized = phase_rotation_identification.*U_identification;
y_phase_normalized = phase_rotation_identification.*y(identification_indices);
U_phase_normalized_selected = U_phase_normalized(:, support_complex);

% Build and structurally reduce the explicit PN-IQ features.
[features_iq_raw, iq_details] = buildPhaseNormalizedIQRegressors( ...
    x, identification_indices, gmp_manager, support_complex);
[features_iq, iq_reduction] = removeStructurallyZeroQFeatures( ...
    features_iq_raw, iq_details.featureMetadata, 1e-12);
builder_error = norm(iq_details.UPhaseNormalized - ...
    U_phase_normalized_selected, 'fro') / ...
    max(1, norm(U_phase_normalized_selected, 'fro'));
support_descriptors = descriptors(support_complex);
zero_q_count = nnz([support_descriptors.QColumnStructurallyZero]);
expected_feature_count = 2*numel(support_complex) - zero_q_count;
assert(isequal(size(features_iq_raw), ...
    [numel(identification_indices), 2*numel(support_complex)]));
assert(isreal(features_iq_raw) && all(isfinite(features_iq_raw), 'all'));
assert(builder_error <= 1e-12);
assert(iq_details.maxCanonicalIError <= 1e-12);
assert(iq_details.maxCanonicalQError <= 1e-12);
assert(iq_reduction.structurallyZeroRemoved == zero_q_count);
assert(iq_reduction.structuralDuplicatesRemoved == 0);
assert(iq_reduction.structuralOppositesRemoved == 0);
assert(size(features_iq, 2) == expected_feature_count);

complex_grid = fitComplexGMPGrid(U_identification, ...
    y(identification_indices), support_complex, [0, 1e-4]);
assert(isequal(complex_grid.support, support_complex));
assert(isequal(size(complex_grid.coefficients), ...
    [numel(support_complex), 2]));
assert(all(isfinite(complex_grid.coefficients), 'all'));

% Check r(n)x(n)=|x(n)| and Q=0 for the current carrier column.
current_linear_index = find(arrayfun(@(item) ...
    item.canonicalGMP && item.carrierLag == 0 && ...
    isempty(item.envelopeLags), descriptors), 1);
probe_rows = full_signal_indices(1:128);
[lag_zero_features, lag_zero_details] = buildPhaseNormalizedIQRegressors( ...
    x, probe_rows, gmp_manager, current_linear_index);
probe_scale = max(1, norm(abs(x(probe_rows)), Inf));
assert(norm(lag_zero_features(:, 1) - abs(x(probe_rows)), Inf) ...
    <= 1e-12*probe_scale);
assert(norm(lag_zero_features(:, 2), Inf) <= 1e-12*probe_scale);
assert(lag_zero_details.maxCarrierI0Error <= 1e-12);
assert(lag_zero_details.maxCarrierQ0 <= 1e-12);

% Compare the complex model with its exact coupled-real representation.
exact_cfg = struct('supportMode', 'reuse_complex_support', ...
    'supportComplex', support_complex, 'normUComplex', complex_fit.normU, ...
    'lambda', 1e-3);
exact_fit = fitPhaseNormGMPReal( ...
    U_phase_normalized, y_phase_normalized, exact_cfg);
yhat_complex = GMP_blockPredict( ...
    x, gmp_manager, support_complex, complex_fit.h_ridge_1e3, ...
    smoke_cfg.blockSize, full_signal_indices);
[yhat_exact, yhat_exact_phase] = predictPhaseNormGMPReal( ...
    x, full_signal_indices, gmp_manager, exact_fit, smoke_cfg.blockSize);
exact_relative_error = norm(yhat_exact - yhat_complex) / ...
    max(norm(yhat_complex), realmin);
assert(isequal(exact_fit.supportComplex, support_complex));
assert(numel(exact_fit.hReal) == 2*numel(support_complex));
assert(isequal(size(yhat_exact_phase), [numel(full_signal_indices), 1]));
assert(exact_relative_error < 1e-9);
assert(abs(nmseComplexDb(y(full_signal_indices), yhat_exact) - ...
    nmse_db(y(full_signal_indices), yhat_exact)) < 1e-10);

% Fit and predict the retained independent PN-IQ model.
independent_fit = fitIndependentIQGMP( ...
    features_iq, y_phase_normalized, 1e-3, ...
    iq_reduction, "Smoke independent PN-IQ-GMP");
[yhat_independent, yhat_independent_phase] = predictIndependentIQGMP( ...
    x, full_signal_indices, gmp_manager, support_complex, ...
    independent_fit, smoke_cfg.blockSize);
assert(independent_fit.effectiveRealFeatures == expected_feature_count);
assert(independent_fit.realParameters == 2*expected_feature_count);
assert(isequal(size(yhat_independent_phase), ...
    [numel(full_signal_indices), 1]));

% Verify the independent I/Q control without phase normalization.
conjugate_index = find(arrayfun(@(item) ...
    item.auxiliaryType == "conjugate", descriptors), 1);
envelope_index = find(arrayfun(@(item) ...
    item.auxiliaryType == "envelope", descriptors), 1);
no_pn_support = [conjugate_index; current_linear_index; envelope_index];
[no_pn_raw, no_pn_details] = buildUnnormalizedIQRegressors( ...
    x, identification_indices, gmp_manager, no_pn_support);
[no_pn_features, no_pn_reduction] = reduceStructuralFeatures( ...
    no_pn_raw, no_pn_details.featureMetadata, 1e-12);
no_pn_fit = fitIndependentIQGMP( ...
    no_pn_features, y(identification_indices), 1e-3, ...
    no_pn_reduction, "Smoke independent I/Q without PN");
no_pn_metrics = evaluateIndependentIQFits( ...
    x, y, full_signal_indices, gmp_manager, no_pn_support, ...
    {no_pn_fit}, "no_phase_normalization", smoke_cfg.blockSize);
assert(size(no_pn_raw, 2) == 6 && size(no_pn_features, 2) == 3);
assert(no_pn_reduction.structurallyZeroRemoved == 1);
assert(no_pn_reduction.structuralDuplicatesRemoved == 1);
assert(no_pn_reduction.structuralOppositesRemoved == 1);
assert(no_pn_fit.realParameters == 6);
assert(all(isfinite(no_pn_metrics.NMSE_dB)));

% Scale the shared-feature 200-parameter control to the fixture.
group_feature_count = min(5, size(features_iq, 2));
group_selection = selectSharedIQFeatures( ...
    features_iq, y_phase_normalized, group_feature_count);
group_reduction = iq_reduction;
group_reduction.keptIndices = iq_reduction.keptIndices( ...
    group_selection.supportFeatures);
group_reduction.effectiveFeatureCount = group_feature_count;
group_fit = fitIndependentIQGMP( ...
    features_iq(:, group_selection.supportFeatures), y_phase_normalized, ...
    1e-3, group_reduction, "Smoke DOMP-selected PN-IQ");
[yhat_group, yhat_group_phase] = predictIndependentIQGMP( ...
    x, full_signal_indices, gmp_manager, support_complex, ...
    group_fit, smoke_cfg.blockSize);
assert(numel(group_selection.supportFeatures) == group_feature_count);
assert(numel(unique(group_selection.supportFeatures)) == group_feature_count);
assert(group_selection.selectionDomainRows == numel(identification_indices));
assert(group_fit.realParameters == 2*group_feature_count);
assert(isequal(size(yhat_group_phase), ...
    [numel(full_signal_indices), 1]));

% Verify the retained four-row analytical cost table.
[complexity, operation_details] = countModelOperations( ...
    gmp_manager.regPopulation, support_complex, iq_reduction);
operation_columns = complexity{:, { ...
    'TotalRealMultiplicationsPerSample', 'TotalRealAdditionsPerSample', ...
    'Divisions', 'SquareRoots', 'Comparisons'}};
assert(height(complexity) == 4);
assert(all(isfinite(operation_columns), 'all'));
assert(all(operation_columns >= 0, 'all'));
assert(all(operation_columns == floor(operation_columns), 'all'));
assert(operation_details.activeComplexRegressors == numel(support_complex));

% Verify the fixed FLOP convention and the expanded/core equivalence.
flop_convention = getFLOPConvention();
[complexity_flops, flop_details] = countModelFLOPs( ...
    complexity, flop_convention);
assert(isequal(flop_convention.FLOPs, [1; 1; 2; 2; 6; 8; 1]));
assert(height(complexity_flops) == 4);
assert(isequal(complexity_flops.CoreFLOPsPerSample, ...
    complexity.TotalRealMultiplicationsPerSample + ...
    complexity.TotalRealAdditionsPerSample));
assert(all(isnan(complexity_flops.EstimatedTotalFLOPsPerSample)));
assert(~flop_details.specialOperationsWeighted);

% Verify the PNNN feature dimension, N12 parameterization, and update budget.
[pnnn_features, pnnn_targets, pnnn_rotation] = buildPhaseNormDataset( ...
    x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
pnnn_features = pnnn_features.';
pnnn_targets = pnnn_targets.';
assert(isequal(size(pnnn_features), [n_samples, 84]));
assert(isequal(size(pnnn_targets), [n_samples, 2]));
assert(numel(pnnn_rotation) == n_samples);
parameter_choice = chooseParameterMatchedHiddenUnits(84, 358);
assert(parameter_choice.hiddenNeurons == 4);
assert(parameter_choice.realParameters == 350);
assert(parameter_choice.absoluteDifference == -8);
n12_parameters = countPNNNParameters(84, 12);
assert(n12_parameters.realParameters == 12*(84+3)+2);
assert(n12_parameters.realParameters == 1046);
pnnn_flops = countPNNNFLOPs( ...
    "PNNN H4 dense", 84, 4, cfg.pnnn.M, cfg.pnnn.orders);
assert(pnnn_flops.NumRealWeights == 344);
assert(pnnn_flops.NumRealBiases == 6);
assert(pnnn_flops.CoreFLOPsPerSample == 876);
assert(pnnn_flops.FLOPsPerSample == 876);
assert(pnnn_flops.DenseExecutionCoreFLOPsPerSample == 876);
assert(pnnn_flops.IdealSparseCoreFLOPsPerSample == 876);
assert(pnnn_flops.NumELUPerSample == 4);
assert(pnnn_flops.NumSqrtPerSample == 14);
assert(pnnn_flops.NumRealDivisionsPerSample == 2);

production_budget = scalePNNNTrainingBudget(491520, 16711, ...
    cfg.training, cfg.pruning);
assert(production_budget.historicalTrainingSamples == 344064);
assert(production_budget.historicalDenseIterationsPerEpoch == 336);
assert(production_budget.currentDenseIterationsPerEpoch == 16);
assert(production_budget.historicalDenseUpdates == 50400);
assert(production_budget.denseMaxEpochs == 3150);
assert(production_budget.denseLearnRateDropPeriod == 105);
assert(production_budget.denseValidationPatience == 1050);
assert(production_budget.currentFineTuneIterationsPerEpoch == 17);
assert(production_budget.fineTuneEpochs == 396);

% Exercise the exact dynamic pruning path on a one-epoch neural fixture.
smoke_training = cfg.training;
smoke_training.maxEpochs = 1;
smoke_training.miniBatchSize = 8;
smoke_training.learnRateDropPeriod = 1;
smoke_training.validationPatience = 1;
smoke_pruning = cfg.pruning;
smoke_pruning.targetActiveTrainableParams = ...
    independent_fit.realParameters;
smoke_pruning.fineTuneEpochs = 2;
smoke_pruning.fineTuneLearnRateDropPeriod = 1;
smoke_runtime_cfg = struct('training', smoke_training, ...
    'pruning', smoke_pruning);

n12_dense_smoke = fitFairPNNNDenseValidation( ...
    pnnn_features, pnnn_targets, pnnn_rotation, y, split.internalTrainIndices, ...
    split.internalValidationIndices, split.identificationIndices, 12, 42, smoke_training);
assert(n12_dense_smoke.parameterCount.realParameters == 1046);
assert(~n12_dense_smoke.selectionContract.testRowsUsed);
assert(n12_dense_smoke.selectionContract.trainingRows == ...
    numel(split.internalTrainIndices));
assert(n12_dense_smoke.selectionContract.validationRows == ...
    numel(split.internalValidationIndices));

dynamic_target = independent_fit.realParameters;
n12_sparse_smoke = pruneAndFineTuneFairPNNN( ...
    n12_dense_smoke, pnnn_features, pnnn_targets, pnnn_rotation, y, ...
    split.internalTrainIndices, split.internalValidationIndices, ...
    split.identificationIndices, dynamic_target, 42, smoke_runtime_cfg);
assert(n12_sparse_smoke.targetActiveParams == dynamic_target);
assert(n12_sparse_smoke.actualActiveParams == dynamic_target);
assert(n12_sparse_smoke.activeBiases == 14);
assert(n12_sparse_smoke.activeWeights == dynamic_target - 14);
assert(n12_sparse_smoke.pruningStats.prunedBiasParams == 0);
assert(n12_sparse_smoke.maskIntegrityAfterPruning.ok);
assert(n12_sparse_smoke.maskIntegrityAfterFineTune.ok);
assert(n12_sparse_smoke.maskIntegrityAfterFineTune.violationCount == 0);
assert(~n12_sparse_smoke.selectionContract.testRowsUsed);
assert(n12_sparse_smoke.fineTuneUpdates == ...
    smoke_pruning.fineTuneEpochs * ...
    ceil(numel(split.internalTrainIndices)/smoke_training.miniBatchSize));

learnables = n12_sparse_smoke.network.Learnables;
for row = 1:height(learnables)
    if lower(string(learnables.Parameter(row))) == "bias"
        assert(all(n12_sparse_smoke.pruningState.masks{row}, 'all'));
    end
end

% Verify the fixed final refit uses every identification row and no holdout.
n12_dense_final = refitFairPNNNDense(pnnn_features, pnnn_targets, ...
    pnnn_rotation, y, identification_indices, full_signal_indices, ...
    12, 42, 1, smoke_training);
assert(n12_dense_final.finalFitSamples == numel(identification_indices));
assert(n12_dense_final.normalizationSamples == ...
    numel(identification_indices));
assert(n12_dense_final.fullSignalSamples == n_samples);
assert(~n12_dense_final.finalContract.fullSignalUsedForTraining);
assert(numel(n12_dense_final.fullSignalPrediction) == n_samples);

n12_sparse_final = refitFairPNNNSparse(n12_dense_final, ...
    pnnn_features, pnnn_targets, pnnn_rotation, y, ...
    identification_indices, full_signal_indices, dynamic_target, ...
    42, 1, smoke_runtime_cfg);
assert(n12_sparse_final.actualActiveParams == dynamic_target);
assert(n12_sparse_final.finalFitSamples == numel(identification_indices));
assert(n12_sparse_final.normalizationSamples == ...
    numel(identification_indices));
assert(n12_sparse_final.fullSignalSamples == n_samples);
assert(~n12_sparse_final.finalContract.fullSignalUsedForTraining);
assert(n12_sparse_final.maskIntegrityAfterFineTune.ok);

sparse_flops = countSparsePNNNFLOPs( ...
    "PNNN N12 sparse", 84, 12, cfg.pnnn.M, cfg.pnnn.orders, ...
    n12_sparse_final.activeWeights, n12_sparse_final.activeBiases);
assert(sparse_flops.DenseExecutionCoreFLOPsPerSample == 2252);
assert(sparse_flops.IdealSparseCoreFLOPsPerSample < ...
    sparse_flops.DenseExecutionCoreFLOPsPerSample);
assert(sparse_flops.FLOPsPerSample == ...
    sparse_flops.IdealSparseCoreFLOPsPerSample);
assert(sparse_flops.DenseMatrixFLOPsPerSample == 2252);
assert(sparse_flops.SparseZeroWeightsSkipped);
assert(sparse_flops.IdealSparseCostRequiresSparseKernel);
assert(sparse_flops.NumRealParameters == dynamic_target);

% Exercise the complete six-model linear study under the new protocol.
study_cfg = cfg;
study_cfg.gmp = smoke_cfg;
study_cfg.gmp.maxPopulation = 6;
study_cfg.gmp.dompOptions.residualTolerance = 0;
study_cfg.gmp.dompOptions.correlationTolerance = 0;
study_cfg.lambdaGrid = [0, 1e-4];
study_cfg.reducedRealParameterTarget = 10;
y_study = y + 1e-3*(randn(n_samples, 1) + 1j*randn(n_samples, 1));
linear_study = runPNGMPDOMPStudy(x, y_study, split, study_cfg);
assert(height(linear_study.comparisonResults) == 6);
assert(all(isfinite(linear_study.comparisonResults.IdentificationNMSEdB)));
assert(all(isfinite(linear_study.comparisonResults.FullSignalNMSEdB)));
assert(all(linear_study.comparisonResults.FLOPsPerSample > 0));
assert(numel(linear_study.supports.supportDOMP) == 6);
assert(numel(linear_study.fullSignalPredictions.complexGMP) == n_samples);
assert(linear_study.equivalenceRelativeError < 1e-9);
independent_row = linear_study.comparisonResults.Model == ...
    "Independent PN-IQ full";
matched_row = linear_study.comparisonResults.Model == ...
    "Complex GMP DOMP parameter-matched";
assert(linear_study.comparisonResults.NumRealParameters(matched_row) == ...
    linear_study.comparisonResults.NumRealParameters(independent_row));

% Every retained coefficient and prediction path must remain finite.
assert(all(isfinite([features_iq(:); no_pn_features(:)])));
assert(all(isfinite([exact_fit.hReal; independent_fit.coefficientsI; ...
    independent_fit.coefficientsQ; no_pn_fit.coefficientsI; ...
    no_pn_fit.coefficientsQ; group_fit.coefficientsI; ...
    group_fit.coefficientsQ])));
assert(all(isfinite([yhat_complex; yhat_exact; yhat_independent; yhat_group])));

fprintf('\nFULL-SIGNAL DOMP COMPARISON SMOKE TEST: PASS\n');
fprintf(['Identification: %d | full signal: %d | population: %d | ' ...
    'support: %d\n'], numel(identification_indices), ...
    numel(full_signal_indices), ...
    n_population, numel(support_complex));
fprintf('Coupled relative prediction error: %.6e\n', exact_relative_error);
fprintf('Independent PN-IQ features/parameters: %d / %d\n', ...
    independent_fit.effectiveRealFeatures, independent_fit.realParameters);
