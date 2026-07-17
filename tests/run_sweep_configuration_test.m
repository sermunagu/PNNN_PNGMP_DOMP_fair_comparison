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
assert(~isfield(cfg.sweep, 'includeHistoricalPNIQReference'));
assert(~isfield(cfg.sweep, 'historicalReferenceParameters'));
assert(isequal(cfg.fixedRidgeLambdas, [1e-3 1e-4 1e-5]));
assert(~isfield(cfg, 'fixedRidgeLambda'));
assert(~isfield(cfg, 'historicalDisjointResultDirectory'));
assert(~isfield(cfg.gmp, 'selectionMethod'));
assert(~isfield(cfg.pnnn, 'denseControlHiddenNeurons'));
assert(~isfield(cfg.pnnn, 'parameterMatchedTargetModel'));
assert(~isfield(cfg.pnnn, 'trainHistoricalN25'));
assert(~isfield(cfg.sweep, 'warmStart'));
assert(~isfield(cfg, 'report'));

fprintf('SWEEP CONFIGURATION TEST: PASS\n');
