% Test the public three-family sweep presentation on synthetic tables.
% The fixture verifies exact matching and treats every requested target equally.
% No model fitting, pruning, measurement loading, or checkpoint reuse is required.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

targets = [300; 344; 380];
Model = [repmat("Complex GMP DOMP sweep", 3, 1); ...
    repmat("Independent PN-IQ PN-DOMP sweep", 3, 1); ...
    repmat("Sparse PNNN N12", 3, 1)];
SweepRole = repmat("Sweep point", 9, 1);
TargetRealParameters = repmat(targets, 3, 1);
ActualRealParameters = TargetRealParameters;
FullSignalNMSEdB = [-31.0; -31.5; -31.8; ...
    -31.2; -31.7; -32.0; -28.0; -28.5; -29.0];
FLOPsPerSample = [900; 1000; 1100; 700; 760; 820; 600; 650; 700];
ActiveWeights = ActualRealParameters;
ActiveBiases = zeros(9, 1);
pnnnRows = Model == "Sparse PNNN N12";
ActiveBiases(pnnnRows) = 14;
ActiveWeights(pnnnRows) = ...
    ActualRealParameters(pnnnRows) - ActiveBiases(pnnnRows);
fixture = table(Model, SweepRole, TargetRealParameters, ...
    ActualRealParameters, FullSignalNMSEdB, FLOPsPerSample, ...
    ActiveWeights, ActiveBiases);

linear = struct('complexTable', fixture(1:3, :), ...
    'pnTable', fixture(4:6, :));
outputDirectory = tempname;
mkdir(outputDirectory);
outputCleanup = onCleanup(@() rmdir(outputDirectory, 's'));

results = writeSweepPresentationOutputs( ...
    linear, fixture(7:9, :), targets, outputDirectory);
assert(height(results) == 3*numel(targets));
assert(~any(results.SweepRole == "Historical reference"));
assert(~any(contains(results.Model, "Historical")));
for index = 1:numel(targets)
    rows = results.TargetRealParameters == targets(index);
    assert(nnz(rows) == 3);
    assert(all(results.ActualRealParameters(rows) == targets(index)));
end
pnnnRows = results.Model == "Sparse PNNN N12";
assert(all(results.ActiveWeights(pnnnRows) + ...
    results.ActiveBiases(pnnnRows) == ...
    results.ActualRealParameters(pnnnRows)));
manualTargetRows = results.TargetRealParameters == 344;
assert(nnz(manualTargetRows) == 3);
assert(all(results.ActualRealParameters(manualTargetRows) == 344));

csv = readtable(fullfile(outputDirectory, 'complexity_sweep.csv'), ...
    'TextType', 'string');
assert(height(csv) == 3*numel(targets));
assert(~any(csv.SweepRole == "Historical reference"));
assert(~any(contains(csv.Model, "Historical")));
assert(isfile(fullfile(outputDirectory, ...
    'comparison_nmse_parameters_sweep.png')));
assert(isfile(fullfile(outputDirectory, ...
    'comparison_nmse_flops_sweep.png')));
clear outputCleanup;

fprintf('SWEEP PRESENTATION TEST: PASS\n');
