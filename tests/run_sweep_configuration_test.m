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
assert(cfg.sweep.schemaVersion == 3);
assert(cfg.selection.nmseToleranceDb == 0.20);
assert(isequal(cfg.selection.sensitivityTolerancesDb, ...
    [0.10 0.15 0.20 0.25]));
assert(cfg.selection.criterionName == ...
    "near-optimal minimum-complexity criterion");

fprintf('SWEEP CONFIGURATION TEST: PASS\n');
