function [modelFile, comparison] = runPNNNPhaseExperiment(sweepDirectory, P)
% Explicit fresh PNNN-only experiment; keep baseline and both linear models intact.
% This runs the FULL historical dense-selection/refit protocol. Not a quick diagnostic.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
data = load(fullfile(sweepDirectory, 'experiment.mat'));
saved = load(fullfile(sweepDirectory, sprintf('models_P%04d.mat', P)), 'bundle');
baseline = saved.bundle;
cfg = data.cfg;
cfg.pnnn.zeroCurrentImaginaryTap = true;
cfg.sweep.parameterGrid = P;
resultsRoot = fullfile(root, 'results', 'pnnn_phase_experiment');
if ~isfolder(resultsRoot), mkdir(resultsRoot); end
directory = tempname(resultsRoot);
mkdir(directory);
fprintf('New zero-current-imaginary-tap PNNN experiment: %s\n', directory);
baselineFile = fullfile(sweepDirectory, sprintf('models_P%04d.mat', P));
save(fullfile(directory, 'experiment.mat'), 'cfg', 'baselineFile', '-v7.3');
timer = tic;
[denseSource, features, targets, rotation] = trainDensePNNN( ...
    data.x, data.y, data.split, cfg, cfg.reducedRealParameterTarget);
save(fullfile(directory, 'dense_source.mat'), 'denseSource', '-v7.3');
denseSource.runtimeConfig.training.verbose = cfg.training.verbose;
point = fitSparsePNNN(denseSource, P, features, targets, ...
    rotation, data.y, data.split, cfg);
elapsedSeconds = toc(timer);
evaluation = setdiff(data.split.fullSignalIndices, data.split.identificationIndices);
point.row.EvaluationNMSEdB = nmseComplexDb(data.y(evaluation), point.fullSignalPrediction(evaluation));
save(fullfile(directory, 'pnnn.mat'), 'point', 'elapsedSeconds', '-v7.3');
bundle = baseline;
bundle.models{3} = point.model;
bundle.config = cfg;
neuralRow = bundle.metrics.Model == cfg.names.pnnn;
if ~ismember('EvaluationNMSEdB', bundle.metrics.Properties.VariableNames)
    bundle.metrics.EvaluationNMSEdB(:) = NaN;
end
comparison = [bundle.metrics(neuralRow,:); point.row];
comparison.Variant = ["Historical baseline"; "Structural zero (fresh training)"];
bundle.metrics(neuralRow,:) = point.row;
bundle.phaseExperiment = struct('baselineFile', baselineFile, ...
    'variant', "zero_current_imaginary_tap_v1", 'elapsedSeconds', elapsedSeconds, ...
    'denseEpochs', denseSource.bestDenseEpoch, 'fineTuneEpochs', denseSource.fineTuneEpochs);
modelFile = fullfile(directory, 'models.mat');
save(modelFile, 'bundle', '-v7.3');
writetable(comparison, fullfile(directory, 'comparison.csv'));
disp(comparison);
fprintf('Models saved: %s\n', modelFile);
end
