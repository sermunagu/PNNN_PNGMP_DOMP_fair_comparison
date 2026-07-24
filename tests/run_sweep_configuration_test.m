% Test the public sweep grid and runner navigation text.
% The fixture reads configuration and source without fitting any model.
% A manually requested 344-point remains an ordinary presentation target.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'config'));

cfg = getFairDOMPComparisonConfig(projectRoot);
assert(isequal(cfg.sweep.parameterGrid, 20:10:500));
assert(~ismember(344, cfg.sweep.parameterGrid));
assert(isequal(cfg.fixedRidgeLambdas, [1e-3 1e-4 1e-5]));
assert(cfg.gmp.dompOptions.columnTolerance == 1e-12);
assert(isequal(cfg.pnnn.orders, [1 3 5 7]));
assert(cfg.pnnn.sparseBaseHiddenNeurons == 12);
assert(cfg.pnnn.nnSeed == 42);
assert(cfg.sweep.schemaVersion == 8);
assert(cfg.names.complexGMPDOMP == "Complex GMP-DOMP");
assert(cfg.names.pniqGMP == "PN-IQ-GMP");
assert(cfg.names.pnnn == "PNNN");
assert(cfg.names.measuredOutput == "Measured output");
assert(cfg.paper.validationNMSELabel == "Validation NMSE (dB)");
assert(cfg.sweep.coefficientRangeDefinition == ...
    "unit_peak_output_per_column_peak_regressors_v4");
assert(cfg.sweep.linearIdentificationScope == ...
    "complete identification subset");
assert(cfg.sweep.linearPrincipalLambda == 0);
assert(cfg.sweep.linearLambdaSelection == "none");
assert(cfg.sweep.fixedRidgeSupportPolicy == ...
    "reuse principal identification DOMP path");
assert(cfg.selection.stabilizationWindowParameters == 100);
assert(cfg.selection.stabilizationToleranceDb == 0.20);
assert(isequal(cfg.selection.sensitivityWindowsParameters, [80 100 120]));
assert(isequal(cfg.selection.sensitivityTolerancesDb, [0.15 0.20 0.25]));
assert(~isfield(cfg.selection, 'criterionName'));

runnerSource = fileread(fullfile(projectRoot, 'run_parameter_sweep.m'));
assert(contains(runnerSource, 'hasCurrentLinearProtocol'));
assert(contains(runnerSource, 'hasCurrentFixedProtocol'));
assert(~contains(runnerSource, '''lambdaGrid'', cfg.lambdaGrid'));
assert(contains(runnerSource, '''names'', cfg.names'));

assert(isfile(fullfile(projectRoot, 'toolbox', 'sweep', 'fit_pniq_gmp.m')));
assert(~isfile(fullfile(projectRoot, 'toolbox', 'sweep', ...
    ['fit_independent_pniq_' 'domp.m'])));

pnnnSource = fileread(fullfile(projectRoot, 'toolbox', 'pnnn', ...
    'prepare_pnnn_dense_source.m'));
assert(contains(pnnnSource, 'split.internalTrainIndices'));
assert(contains(pnnnSource, 'split.internalValidationIndices'));

fprintf('SWEEP CONFIGURATION TEST: PASS\n');
