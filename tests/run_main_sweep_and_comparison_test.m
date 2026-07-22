% Test automatic selection and explicit manual override with shadowed I/O.
% No measurement, sweep artifact, DOMP path, or neural model is used.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
fixtureDirectory = tempname;
mkdir(fixtureDirectory);
originalDirectory = pwd;
warningState = warning('off', 'MATLAB:dispatcher:nameConflict');
cleanup = onCleanup(@() restoreFixture( ...
    originalDirectory, fixtureDirectory, warningState));

writelines([ ...
    "function sweep = run_parameter_sweep(targets)"; ...
    "global mockSweepTargets"; ...
    "mockSweepTargets = targets;"; ...
    "selection = struct('selectedParameters', 360, " + ...
        "'summarySentence', 'synthetic quantified selection');"; ...
    "sweep = struct('fixture', true, 'selection', selection);"; ...
    "end"], fullfile(fixtureDirectory, 'run_parameter_sweep.m'));
writelines([ ...
    "function results = run_selected_comparison(target, sweep)"; ...
    "global mockSelectedTargets"; ...
    "mockSelectedTargets(end+1) = target;"; ...
    "results = struct('selectedParameters', target, 'sweep', sweep);"; ...
    "end"], fullfile(fixtureDirectory, ...
    'run_selected_comparison.m'));

addpath(projectRoot, '-end');
cd(fixtureDirectory);
rehash;
clear run_parameter_sweep run_selected_comparison ...
    main_sweep_and_comparison;
global mockSweepTargets mockSelectedTargets %#ok<GVMIS>
mockSweepTargets = [];
mockSelectedTargets = [];

automatic = main_sweep_and_comparison();
assert(automatic.selectedParameters == 360);
assert(automatic.selection.selectedParameters == 360);
assert(isequal(mockSweepTargets, 20:10:500));
manual = main_sweep_and_comparison(350);
assert(manual.selectedParameters == 350);
assert(manual.selection.selectedParameters == 360);
assert(isequal(mockSelectedTargets, [360 350]));

selectedSource = fileread(fullfile(projectRoot, ...
    'run_selected_comparison.m'));
assert(contains(selectedSource, 'buildCompleteComparisonTable'));
assert(contains(selectedSource, 'selected_complete_comparison.csv'));
assert(contains(selectedSource, "'completeComparisonTable'"));
assert(contains(selectedSource, 'selected_time_domain_real'));
assert(contains(selectedSource, 'selected_time_domain_imaginary'));
assert(contains(selectedSource, 'timeDomainFigureFiles'));
assert(contains(selectedSource, ...
    'exportSelectedTimeDomainFigures(targetFullSignal, complexPrediction'));
assert(contains(selectedSource, ...
    'pniqPrediction, selectedDirectory, cfg.names, exportOptions'));
assert(~contains(selectedSource, 'bestPNNNParameters'));

clear cleanup;
fprintf('MAIN SWEEP AND COMPARISON TEST: PASS\n');

function restoreFixture(directory, fixtureDirectory, warningState)
cd(directory);
clear run_parameter_sweep run_selected_comparison ...
    main_sweep_and_comparison;
warning(warningState);
if isfolder(fixtureDirectory)
    rmdir(fixtureDirectory, 's');
end
end
