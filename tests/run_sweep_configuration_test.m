% Test the public sweep grid, minimum N12 target, and runner navigation text.
% The fixture stops before measurement loading and performs no model fitting.
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

try
    run_parameter_sweep(14);
    error('run_sweep_configuration_test:ExpectedMinimumError', ...
        'A target containing only the protected biases must be rejected.');
catch exception
    assert(exception.identifier == "run_parameter_sweep:InvalidTargets");
end

source = string(fileread(fullfile(projectRoot, 'run_parameter_sweep.m')));
for forbidden = ["rows344", "pnComparison344", ...
        "historicalReferenceParameters", "includeHistoricalPNIQReference"]
    assert(~contains(source, forbidden));
end
for required = ["Signed parameter-complexity sweep", "[Linear]", ...
        "[PNNN", "[Output]", "Sweep completed"]
    assert(contains(source, required));
end

fprintf('SWEEP CONFIGURATION TEST: PASS\n');
