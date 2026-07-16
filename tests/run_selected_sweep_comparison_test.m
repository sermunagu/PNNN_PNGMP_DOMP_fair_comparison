% Test the final three-family comparison with tiny signed checkpoints.
% The fixture validates exact parameters, signatures, reuse, and predictions.
% It performs no measurement loading, DOMP, fitting, pruning, or training.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
fixtureDirectory = tempname;
mkdir(fixtureDirectory);
fixtureCleanup = onCleanup(@() rmdir(fixtureDirectory, 's'));

target = 340;
Model = ["Complex GMP DOMP sweep"; ...
    "Independent PN-IQ PN-DOMP sweep"; "Sparse PNNN N12"];
SweepRole = repmat("Sweep point", 3, 1);
TargetRealParameters = repmat(target, 3, 1);
ActualRealParameters = TargetRealParameters;
FullSignalNMSEdB = [-31.4; -31.8; -29.1];
FLOPsPerSample = [1400; 820; 790];
ActiveWeights = [target; target; target - 14];
ActiveBiases = [0; 0; 14];
rows = table(Model, SweepRole, TargetRealParameters, ...
    ActualRealParameters, FullSignalNMSEdB, FLOPsPerSample, ...
    ActiveWeights, ActiveBiases);

signature = struct('schemaVersion', 2, 'algorithm', ...
    "fixture-signature", 'digest', "fixture-experiment");
identity = struct('schemaVersion', 2, 'digest', "fixture-sweep");
sweep = struct('results', rows, 'resultDirectory', ...
    string(fixtureDirectory), 'experimentSignature', signature, ...
    'sweepIdentity', identity);

n = 24;
complexPrediction = complex((1:n).', -(1:n).');
pnPrediction = 0.9*complexPrediction;
pnnnPrediction = 0.8*complexPrediction;
supports = struct('complex', {{(1:target/2).'}}, ...
    'pnFeatures', {{(1:target/2).'}}, ...
    'pnComplex', {{(1:target/2).'}});
linearPayload = struct('complexTable', rows(1,:), 'pnTable', rows(2,:), ...
    'supports', supports, 'paths', struct(), ...
    'predictions', struct('complexFull', complexPrediction, ...
    'pnFull', pnPrediction));
denseSignature = struct('digest', "fixture-network");
densePayload = struct('signature', denseSignature, 'denseFit', struct());
pnnnPayload = struct('target', target, 'row', rows(3,:), ...
    'fullSignalPrediction', pnnnPrediction, ...
    'denseSourceSignature', denseSignature);
summaryPayload = struct('results', rows);

files = ["linear_sweep.mat", "sweep_dense_source.mat", ...
    "pnnn_target_0340.mat", "complexity_sweep.mat"];
payloads = {linearPayload, densePayload, pnnnPayload, summaryPayload};
for index = 1:numel(files)
    writeArtifact(fullfile(fixtureDirectory, files(index)), ...
        identity, signature, payloads{index});
end
bytesBefore = cellfun(@(name) readBinaryFile(fullfile( ...
    fixtureDirectory, name)), cellstr(files), 'UniformOutput', false);

results = run_fair_PNNN_vs_PNGMP_DOMP(target, sweep);
assert(results.selectedParameters == target);
assert(height(results.comparisonTable) == 3);
assert(isequal(results.comparisonTable.Model, ...
    ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"; "Sparse PNNN N12"]));
assert(all(results.comparisonTable.ActualRealParameters == target));
assert(results.reusedLinearSweep && results.reusedPNNNPoint);
assert(results.reusedDenseSource);
assert(isequal(results.fullSignalPredictions.complexGMP, complexPrediction));
assert(isequal(results.fullSignalPredictions.pnIQ, pnPrediction));
assert(isequal(results.fullSignalPredictions.sparsePNNNN12, pnnnPrediction));
assert(~any(contains(results.comparisonTable.Model, ...
    ["Historical", "H4", "dense"], 'IgnoreCase', true), 'all'));

expectError(@() run_fair_PNNN_vs_PNGMP_DOMP(), ...
    "run_fair_PNNN_vs_PNGMP_DOMP:MissingTarget");
for invalid = {NaN, Inf, 14, 340.5}
    expectError(@() run_fair_PNNN_vs_PNGMP_DOMP(invalid{1}, sweep), ...
        invalidIdentifier(invalid{1}));
end
for ordinaryMissingTarget = [200 344]
    expectError(@() run_fair_PNNN_vs_PNGMP_DOMP( ...
        ordinaryMissingTarget, sweep), ...
        "run_fair_PNNN_vs_PNGMP_DOMP:MissingSweepPoint");
end
incompatibleSweep = sweep;
incompatibleSweep.experimentSignature.digest = "different-experiment";
expectError(@() run_fair_PNNN_vs_PNGMP_DOMP(target, incompatibleSweep), ...
    "run_fair_PNNN_vs_PNGMP_DOMP:IncompatibleArtifact");

source = string(fileread(fullfile(projectRoot, ...
    'run_fair_PNNN_vs_PNGMP_DOMP.m')));
for forbidden = ["runPNGMPDOMPStudy", "runPNNNComparisonStudy", ...
        "selectDOMPSupport", "fitFair", "344", "200"]
    assert(~contains(source, forbidden));
end
bytesAfter = cellfun(@(name) readBinaryFile(fullfile( ...
    fixtureDirectory, name)), cellstr(files), 'UniformOutput', false);
assert(isequal(bytesAfter, bytesBefore));

clear fixtureCleanup;
fprintf('SELECTED SWEEP COMPARISON TEST: PASS\n');

function writeArtifact(filename, identity, signature, payload)
checkpointArtifact = struct('schemaVersion', 1, ...
    'sweepIdentity', identity, 'experimentSignature', signature, ...
    'payload', payload);
save(char(filename), 'checkpointArtifact');
end

function bytes = readBinaryFile(filename)
file = fopen(filename, 'rb');
assert(file >= 0);
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, Inf, '*uint8');
clear cleanup;
end

function expectError(action, expectedIdentifier)
try
    action();
    error('run_selected_sweep_comparison_test:ExpectedError', ...
        'Expected error %s was not raised.', expectedIdentifier);
catch exception
    assert(string(exception.identifier) == expectedIdentifier);
end
end

function identifier = invalidIdentifier(value)
if isequal(value, 14)
    identifier = "run_fair_PNNN_vs_PNGMP_DOMP:TargetBelowPNNNMinimum";
else
    identifier = "run_fair_PNNN_vs_PNGMP_DOMP:InvalidTarget";
end
end
