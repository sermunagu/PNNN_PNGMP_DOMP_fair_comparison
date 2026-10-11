function [modelFile, bundle] = trainDPDOffline(measurementFile, P, cfg)
% Fit three DPD models from scratch and save them for ONLINE.
% The optional cfg is used by the small synthetic tests.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'toolbox')));

if nargin < 2, P = 340; end
if nargin < 3, cfg = getFairDOMPComparisonConfig(root); end

cfg.measurementFile = measurementFile;
cfg.sweep.parameterGrid = P;

%% Load and prepare the selected capture
S = load(measurementFile);
x = double(S.x(:));
y = double(S.y(:));

if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end

split = buildCommonComparisonSplit(x, y, cfg);

%% Fresh DOMP fits and dense PNNN, followed by pruning and fine-tuning
[linear, models] = run_linear_sweep(x, y, split, cfg);

[denseSource, features, targets, rotation] = prepare_pnnn_dense_source( ...
    x, y, split, cfg, cfg.reducedRealParameterTarget);

denseSource.runtimeConfig.training.verbose = cfg.training.verbose;

point = fit_sparse_pnnn_target(denseSource, P, features, targets, rotation, y, split, cfg);

models{3} = point.model;
metrics = [linear.complexTable; linear.pniqTable; point.row];

%% Save one bundle and the direct predictions in a new folder
resultsRoot = fullfile(root, 'results', 'dpd_offline');

if isfield(cfg, 'dpd'), resultsRoot = cfg.dpd.resultsRoot; end
if ~isfolder(resultsRoot), mkdir(resultsRoot); end

outDir = tempname(resultsRoot);
mkdir(outDir);
bundle = struct('P', P, 'models', {models}, 'config', cfg, ...
    'sourceFile', measurementFile, 'split', split, 'metrics', metrics);

if isfield(S, 'fs'), bundle.fs = S.fs; end

modelFile = fullfile(outDir, 'models.mat');

save(modelFile, 'bundle', '-v7.3');
predictions = [linear.predictions.complexFull, linear.predictions.pniqFull, point.fullSignalPrediction];
save(fullfile(outDir, 'predictions.mat'), 'predictions', '-v7.3');
fprintf('Models saved: %s\n', modelFile);

end
