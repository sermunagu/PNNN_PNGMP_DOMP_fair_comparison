% Test deterministic experiment signatures.
% The fixture changes data and science-critical configuration without training.

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
assert(signature.schemaVersion == 3);

x_modified = x;
x_modified(37) = x_modified(37) + eps(x_modified(37));
assert(buildExperimentSignature(x_modified, y, cfg).digest ~= ...
    signature.digest);
assert(buildExperimentSignature(y, x, cfg).digest ~= signature.digest);

cfg_changed = cfg;
cfg_changed.pnnn.M = cfg.pnnn.M + 1;
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.pnnn.orders = [1 3 5];
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.pnnn.nnSeed = cfg.pnnn.nnSeed + 1;
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.training.learnRateDropFactor = 0.9;
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);
cfg_changed = cfg;
cfg_changed.pruning.fineTuneSeedOffset = ...
    cfg.pruning.fineTuneSeedOffset + 1;
assert(buildExperimentSignature(x, y, cfg_changed).digest ~= ...
    signature.digest);

fprintf('EXPERIMENT SIGNATURE TEST: PASS\n');
