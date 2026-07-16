% Test canonical and supplementary sweep presentation on synthetic tables.
% The full 49-target fixture verifies row counts and nine parameter curves.
% No model fitting, DOMP selection, PNNN training, or measurements are used.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

targets = (20:10:500).';
nTargets = numel(targets);
Model = [repmat("Complex GMP DOMP sweep", nTargets, 1); ...
    repmat("Independent PN-IQ PN-DOMP sweep", nTargets, 1); ...
    repmat("Sparse PNNN N12", nTargets, 1)];
SweepRole = repmat("Sweep point", 3*nTargets, 1);
TargetRealParameters = repmat(targets, 3, 1);
ActualRealParameters = TargetRealParameters;
FullSignalNMSEdB = [-25 - 0.01*targets; ...
    -25.5 - 0.01*targets; -22 - 0.008*targets];
FLOPsPerSample = [700 + targets; 500 + targets; 400 + targets];
ActiveWeights = ActualRealParameters;
ActiveBiases = zeros(3*nTargets, 1);
pnnnMask = Model == "Sparse PNNN N12";
ActiveBiases(pnnnMask) = 14;
ActiveWeights(pnnnMask) = ...
    ActualRealParameters(pnnnMask) - ActiveBiases(pnnnMask);
fixture = table(Model, SweepRole, TargetRealParameters, ...
    ActualRealParameters, FullSignalNMSEdB, FLOPsPerSample, ...
    ActiveWeights, ActiveBiases);

linear = struct('complexTable', fixture(1:nTargets, :), ...
    'pnTable', fixture(nTargets + (1:nTargets), :));
lambdas = [1e-3; 1e-4; 1e-5];
variantTargets = repelem(targets, numel(lambdas));
variantLambdas = repmat(lambdas, nTargets, 1);
fixedModel = [repmat("Complex GMP-DOMP", numel(variantTargets), 1); ...
    repmat("PN-IQ PN-DOMP", numel(variantTargets), 1)];
TargetRealParameters = [variantTargets; variantTargets];
ActualRealParameters = TargetRealParameters;
FixedLambda = [variantLambdas; variantLambdas];
IdentificationNMSEdB = -24 - 0.01*ActualRealParameters;
FullSignalNMSEdB = IdentificationNMSEdB - 0.2;
FLOPsPerSample = [repelem(linear.complexTable.FLOPsPerSample, 3); ...
    repelem(linear.pnTable.FLOPsPerSample, 3)];
fixedLinear.table = table(fixedModel, TargetRealParameters, ...
    ActualRealParameters, FixedLambda, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, 'VariableNames', ...
    {'Model','TargetRealParameters','ActualRealParameters','FixedLambda', ...
    'IdentificationNMSEdB','FullSignalNMSEdB','FLOPsPerSample'});

outputDirectory = tempname;
mkdir(outputDirectory);
outputCleanup = onCleanup(@() rmdir(outputDirectory, 's'));
linearCheckpoint = fullfile(outputDirectory, 'linear_sweep.mat');
checkpointSentinel = uint8(0:31);
save(linearCheckpoint, 'checkpointSentinel');
checkpointBytesBefore = readBinaryFile(linearCheckpoint);

[results, details] = writeSweepPresentationOutputs( ...
    linear, fixedLinear, fixture(2*nTargets + (1:nTargets), :), ...
    targets, outputDirectory);
checkpointBytesAfter = readBinaryFile(linearCheckpoint);
assert(isequal(checkpointBytesAfter, checkpointBytesBefore));
assert(height(results) == 147);
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
assert(details.parameterCurveCount == 9);
assert(details.flopsCurveCount == 3);

csv = readtable(fullfile(outputDirectory, 'complexity_sweep.csv'), ...
    'TextType', 'string');
fixedCSV = readtable(fullfile(outputDirectory, ...
    'fixed_lambda_linear_sweep.csv'), 'TextType', 'string');
assert(height(csv) == 147);
assert(height(fixedCSV) == 294);
assert(numel(unique(fixedCSV.Model)) == 2);
for model = unique(fixedCSV.Model).'
    rows = fixedCSV.Model == model;
    assert(isequal(sort(unique(fixedCSV.FixedLambda(rows))), sort(lambdas)));
    for lambda = lambdas.'
        assert(nnz(rows & fixedCSV.FixedLambda == lambda) == nTargets);
    end
end
assert(isfile(fullfile(outputDirectory, ...
    'comparison_nmse_parameters_sweep.png')));
assert(isfile(fullfile(outputDirectory, ...
    'comparison_nmse_flops_sweep.png')));
clear outputCleanup;

fprintf('SWEEP PRESENTATION TEST: PASS\n');

function bytes = readBinaryFile(filename)
file = fopen(filename, 'rb');
assert(file >= 0);
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, Inf, '*uint8');
clear cleanup;
end
