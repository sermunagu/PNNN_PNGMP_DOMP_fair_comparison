% Test silent loading of optional experiment signatures from legacy MAT files.
% Missing signatures return an empty structure without weakening validation.
% The fixture creates only temporary MAT files and performs no model fitting.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));
fixtureDirectory = tempname;
mkdir(fixtureDirectory);
fixtureCleanup = onCleanup(@() rmdir(fixtureDirectory, 's'));

legacyFile = fullfile(fixtureDirectory, 'legacy.mat');
legacy_value = 1;
save(legacyFile, 'legacy_value');
lastwarn('');
legacy = loadOptionalExperimentSignature(legacyFile);
[warningMessage, warningIdentifier] = lastwarn;
assert(isempty(fieldnames(legacy)));
assert(isempty(warningMessage) && isempty(warningIdentifier));

signedFile = fullfile(fixtureDirectory, 'signed.mat');
experiment_signature = struct('schemaVersion', 2, ...
    'digest', "fixture");
save(signedFile, 'experiment_signature');
signed = loadOptionalExperimentSignature(signedFile);
assert(isequaln(signed.experiment_signature, experiment_signature));
assert(isempty(fieldnames(loadOptionalExperimentSignature( ...
    fullfile(fixtureDirectory, 'missing.mat')))));

try
    loadOptionalExperimentSignature(["first.mat", "second.mat"]);
    error('run_optional_signature_load_test:ExpectedError', ...
        'A string array path must be rejected.');
catch exception
    assert(exception.identifier == ...
        "loadOptionalExperimentSignature:InvalidPath");
end

clear fixtureCleanup;
fprintf('OPTIONAL SIGNATURE LOAD TEST: PASS\n');
