% Verify joint stabilization, minimum-budget selection, sensitivity, and regression.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

budgets = (100:10:150).';
gmpNMSE = [-35.00; -35.20; -35.30; -35.35; -35.37; -35.38];
pnNMSE = [-36.00; -36.40; -36.70; -36.80; -36.85; -36.87];
pnnnNMSE = [-30.00; -31.00; -32.00; -32.15; -32.18; -32.19];
gmpFLOPs = [500; 520; 540; 560; 580; 600];
pnFLOPs = [400; 415; 430; 445; 460; 475];
pnnnFLOPs = [300; 320; 340; 360; 380; 400];
Model = [repmat("Complex GMP DOMP sweep", numel(budgets), 1); ...
    repmat("Independent PN-IQ PN-DOMP sweep", numel(budgets), 1); ...
    repmat("Sparse PNNN N12", numel(budgets), 1)];
ActualRealParameters = repmat(budgets, 3, 1);
FullSignalNMSEdB = [gmpNMSE; pnNMSE; pnnnNMSE];
FLOPsPerSample = [gmpFLOPs; pnFLOPs; pnnnFLOPs];
results = table(Model, ActualRealParameters, FullSignalNMSEdB, ...
    FLOPsPerSample);
config = struct('stabilizationWindowParameters', 20, ...
    'stabilizationToleranceDb', 0.20, ...
    'sensitivityWindowsParameters', [10 20 30], ...
    'sensitivityTolerancesDb', [0.10 0.20]);
selection = selectOperatingPoint(results, config);
assert(selection.selectedParameters == 120);
assert(selection.selectedParameters ~= 340);
assert(selection.criterionName == ...
    "joint stabilization minimum-complexity criterion");
assert(selection.selectedComplexGMPFutureGainDb <= 0.20);
assert(selection.selectedPNIQFutureGainDb <= 0.20);
assert(selection.selectedSparsePNNNFutureGainDb <= 0.20);
assert(selection.selectedWorstFutureGainDb <= 0.20);
diagnostics = selection.diagnosticsTable;
assert(~diagnostics.JointlyStabilized( ...
    diagnostics.ActualRealParameters == 110));
assert(diagnostics.JointlyStabilized( ...
    diagnostics.ActualRealParameters == 120));
assert(diagnostics.Selected( ...
    diagnostics.ActualRealParameters == 120));
assert(~diagnostics.HasFullWindow( ...
    diagnostics.ActualRealParameters == 140));

sensitivity = selection.sensitivityTable;
row = sensitivity.StabilizationWindowParameters == 20 & ...
    sensitivity.StabilizationToleranceDb == 0.20;
assert(nnz(row) == 1);
assert(sensitivity.HasJointlyStabilizedPoint(row));
assert(sensitivity.SelectedParameters(row) == 120);

nonmonotonic = results;
row = nonmonotonic.Model == "Independent PN-IQ PN-DOMP sweep" & ...
    nonmonotonic.ActualRealParameters == 130;
nonmonotonic.FLOPsPerSample(row) = 400;
assertError(@() selectOperatingPoint(nonmonotonic, config), ...
    'selectOperatingPoint:NonmonotonicFLOPs');

historicalFile = fullfile(projectRoot, 'results', 'parameter_sweep', ...
    'sweep_d113e389ab78', 'complexity_sweep.csv');
historical = readtable(historicalFile);
historicalConfig = struct('stabilizationWindowParameters', 100, ...
    'stabilizationToleranceDb', 0.20, ...
    'sensitivityWindowsParameters', [80 100 120], ...
    'sensitivityTolerancesDb', [0.15 0.20 0.25]);
historicalSelection = selectOperatingPoint(historical, historicalConfig);
assert(historicalSelection.selectedParameters == 340);
historicalDiagnostics = historicalSelection.diagnosticsTable;
row330 = historicalDiagnostics.ActualRealParameters == 330;
row340 = historicalDiagnostics.ActualRealParameters == 340;
assert(~historicalDiagnostics.JointlyStabilized(row330));
assert(historicalDiagnostics.PNIQFutureGainDb(row330) > 0.20);
assert(historicalDiagnostics.JointlyStabilized(row340));
assert(all([historicalDiagnostics.ComplexGMPFutureGainDb(row340), ...
    historicalDiagnostics.PNIQFutureGainDb(row340), ...
    historicalDiagnostics.SparsePNNNFutureGainDb(row340)] <= 0.20));

fprintf('OPERATING POINT SELECTION TEST: PASS\n');

function assertError(functionHandle, expectedIdentifier)
try
    functionHandle();
catch exception
    assert(strcmp(exception.identifier, expectedIdentifier));
    return;
end
error('run_operating_point_selection_test:MissingError', ...
    'Expected error %s was not raised.', expectedIdentifier);
end
