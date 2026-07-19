% Verify admissibility, complexity tie-breaking, marginals, and regression.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

budgets = [200; 250; 280; 300; 500];
pnnnNMSE = [-39.70; -39.82; -39.85; -39.80; -40.00];
gmpNMSE = [-39.60; -39.90; -39.80; -39.75; -39.70];
pnnnFLOPs = [600; 650; 700; 750; 1000];
gmpFLOPs = 800 + budgets;
Model = [repmat("Complex GMP DOMP sweep", numel(budgets), 1); ...
    repmat("Sparse PNNN N12", numel(budgets), 1)];
ActualRealParameters = [budgets; budgets];
FullSignalNMSEdB = [gmpNMSE; pnnnNMSE];
FLOPsPerSample = [gmpFLOPs; pnnnFLOPs];
results = table(Model, ActualRealParameters, FullSignalNMSEdB, ...
    FLOPsPerSample);
config = struct('nmseToleranceDb', 0.20, ...
    'sensitivityTolerancesDb', [0.10 0.15 0.20 0.25]);
selection = selectOperatingPoint(results, config);
assert(selection.selectedParameters == 280);
assert(selection.selectedParameters ~= 340);
assert(selection.criterionName == ...
    "near-optimal minimum-complexity criterion");
diagnostics = selection.diagnosticsTable;
assert(~diagnostics.BeatsGMP(diagnostics.ActualRealParameters == 250));
assert(diagnostics.Admissible(diagnostics.ActualRealParameters == 280));
assert(diagnostics.Admissible(diagnostics.ActualRealParameters == 300));
assert(diagnostics.Selected(diagnostics.ActualRealParameters == 280));
assert(diagnostics.MarginalGainDbPer100AdditionalFLOPs( ...
    diagnostics.ActualRealParameters == 300) < 0);
tieResults = results;
tieRow = tieResults.Model == "Sparse PNNN N12" & ...
    tieResults.ActualRealParameters == 300;
tieResults.FLOPsPerSample(tieRow) = 700;
tieSelection = selectOperatingPoint(tieResults, config);
assert(tieSelection.selectedParameters == 280);

historicalFile = fullfile(projectRoot, 'results', 'parameter_sweep', ...
    'sweep_0dd97cdd1cca', 'complexity_sweep.csv');
historical = readtable(historicalFile);
historicalSelection = selectOperatingPoint(historical, config);
assert(historicalSelection.nmseToleranceDb == 0.20);
assert(historicalSelection.selectedParameters == 340);
assert(historicalSelection.criterionName == ...
    "near-optimal minimum-complexity criterion");

fprintf('OPERATING POINT SELECTION TEST: PASS\n');
