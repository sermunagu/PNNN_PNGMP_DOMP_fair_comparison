% Test deterministic experiment signatures and safe PNNN artifact reuse.
% The fixture changes data and reuse-critical configuration without training.

clearvars;
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'toolbox', 'pnnn'));

cfg = getFairDOMPComparisonConfig(project_root);
rng(731, 'twister');
x = complex(randn(128, 1), randn(128, 1));
y = complex(randn(128, 1), randn(128, 1));

signature = buildExperimentSignature(x, y, cfg);
signature_repeat = buildExperimentSignature(x, y, cfg);
assert(signature.digest == signature_repeat.digest);

x_modified = x;
x_modified(37) = x_modified(37) + eps(x_modified(37));
assert(buildExperimentSignature(x_modified, y, cfg).digest ~= ...
    signature.digest);
assert(buildExperimentSignature(y, x, cfg).digest ~= signature.digest);

cfg_changed = cfg;
cfg_changed.pnnn.actType = 'tanh';
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.pnnn.M = cfg.pnnn.M + 1;
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.pnnn.orders = [1 3 5];
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);

split.internalTrainIndices = (1:70).';
split.internalValidationIndices = (71:90).';
split.identificationIndices = (1:90).';
savedSplit.internal_train_indices = split.internalTrainIndices;
savedSplit.internal_validation_indices = split.internalValidationIndices;
savedSplit.identification_indices = split.identificationIndices;
savedConfig.experiment_signature = signature;
[matches, reason] = isReusablePNNNSelection( ...
    savedSplit, savedConfig, split, signature);
assert(matches && reason == "compatible");

savedConfig.experiment_signature = ...
    buildExperimentSignature(x_modified, y, cfg);
[matches, reason] = isReusablePNNNSelection( ...
    savedSplit, savedConfig, split, signature);
assert(~matches && reason == "experiment signature mismatch");

savedConfig = struct();
[matches, reason] = isReusablePNNNSelection( ...
    savedSplit, savedConfig, split, signature);
assert(~matches && reason == "missing experiment signature");

fprintf('EXPERIMENT SIGNATURE TEST: PASS\n');
