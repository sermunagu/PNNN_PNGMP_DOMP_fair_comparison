function [bundle, dpd] = testDPD(P, mode)
% Offline ILC imitation. Default: load only; train/recover require an explicit call.
% testDPD(340, "recover") reuses legacy supports and the saved dense PNNN.
if nargin < 1, P = 340; end
if nargin < 2, mode = "reuse"; end
validateattributes(P, {'numeric'}, {'scalar','integer','positive','even'});
mode = validatestring(mode, {'reuse','train','recover'});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'toolbox')));
cfg = getFairDOMPComparisonConfig(root);
cfg.sweep.parameterGrid = P;

%% 1. Desired signal u -> ILC predistorted signal (never PA forward data)
xyFile = fullfile(root, 'measurements', 'experiment20260429T134032_xy.mat');
source = load(xyFile, 'x', 'y', 'fs');
u = source.x(:);
xILC = source.y(:);
assert(numel(u) == numel(xILC) && all(isfinite([u; xILC])));
assert(string(cfg.mappingMode) == "xy_forward", 'DPD requires u -> xILC.');
[x, y] = selectXYByMapping(u, xILC, cfg.mappingMode);
trainingMeans = struct('input', mean(x), 'target', mean(y));
if cfg.pnnn.removeDC
    x = x - trainingMeans.input;
    y = y - trainingMeans.target;
end
split = buildCommonComparisonSplit(x, y, cfg);
signature = buildExperimentSignature(u, xILC, cfg);
tag = char(signature.digest);
outDir = fullfile(root, 'results', 'dpd_ilc', sprintf('P%d_%s', P, tag(1:12)));
cfg.measurementFile = xyFile;
modelFile = fullfile(outDir, sprintf('models_P%d.mat', P));

%% 2. Load final models without identification, or explicitly build missing ones
if isfile(modelFile)
    saved = load(modelFile, 'bundle');
    bundle = saved.bundle;
    assert(isequaln(bundle.config, cfg) && bundle.signature.digest == signature.digest, ...
        'Saved DPD configuration differs. Do not overwrite this experiment.');
else
    if strcmp(mode, 'reuse')
        error('testDPD:MissingModels', ...
            ['Reusable model parameters are missing. Predictions alone are not models. ' ...
            'After authorization use testDPD(%d, "recover") for legacy checkpoints.'], P);
    end
    linearFile = fullfile(outDir, sprintf('linear_P%d.mat', P));
    pnnnFile = fullfile(outDir, sprintf('pnnn_P%d.mat', P));
    denseFile = fullfile(outDir, 'pnnn_dense_source.mat');
    recovering = strcmp(mode, 'recover');
    if recovering
        assert(isfile(linearFile) && isfile(pnnnFile) && isfile(denseFile), ...
            'Recovery needs both legacy checkpoints and the existing dense PNNN.');
        oldLinear = load(linearFile, 'linear');
        oldPNNN = load(pnnnFile, 'point');
    else
        assert(~isfile(linearFile) && ~isfile(pnnnFile), ...
            'Existing predictions found: use recover, not a new identification.');
    end
    if ~isfolder(outDir), mkdir(outDir); end

    %% 3. Linear fits: retain coefficients; recovery keeps the exact DOMP paths
    linearModelsFile = fullfile(outDir, sprintf('linear_models_P%d.mat', P));
    if isfile(linearModelsFile)
        saved = load(linearModelsFile);
        assert(isequaln(saved.config, cfg), 'Linear recovery configuration differs.');
        linear = saved.linear;
        models = saved.models;
    else
        if recovering
            [linear, models] = run_linear_sweep(x, y, split, cfg, oldLinear.linear.paths);
        else
            [linear, models] = run_linear_sweep(x, y, split, cfg);
        end
        config = cfg;
        save(linearModelsFile, 'linear', 'models', 'config', '-v7.3');
    end
    if recovering
        assert(isequal(linear.paths, oldLinear.linear.paths));
        assert(isequal(linear.pniqPathMap, oldLinear.linear.pniqPathMap));
        checkCosts(linear.complexTable, oldLinear.linear.complexTable);
        checkCosts(linear.pniqTable, oldLinear.linear.pniqTable);
        checkPrediction(linear.predictions.complexFull, oldLinear.linear.predictions.complexFull, 1e-10);
        checkPrediction(linear.predictions.pniqFull, oldLinear.linear.predictions.pniqFull, 1e-10);
    end

    %% 4. Sparse PNNN: keep the FINAL network, masks and dense normalization
    sparseModelFile = fullfile(outDir, sprintf('pnnn_model_P%d.mat', P));
    if isfile(sparseModelFile)
        saved = load(sparseModelFile);
        assert(isequaln(saved.config, cfg), 'PNNN recovery configuration differs.');
        point = saved.point;
    elseif recovering && isfield(oldPNNN.point, 'model')
        point = oldPNNN.point;
    else
        if isfile(denseFile)
            saved = load(denseFile, 'denseSource');
            denseSource = saved.denseSource;
            [features, targets, rotation] = buildPhaseNormDataset( ...
                x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
            features = features.';
            targets = targets.';
            rotation = rotation(:);
        else
            % Only explicit train mode can reach dense training.
            assert(~recovering, 'Recovery must never repeat dense training.');
            [denseSource, features, targets, rotation] = ...
                prepare_pnnn_dense_source(x, y, split, cfg, cfg.reducedRealParameterTarget);
            save(denseFile, 'denseSource', '-v7.3');
        end
        denseSource.runtimeConfig.training.verbose = cfg.training.verbose;
        point = fit_sparse_pnnn_target(denseSource, P, ...
            features, targets, rotation, y, split, cfg);
        config = cfg;
        save(sparseModelFile, 'point', 'config', '-v7.3');
    end
    if recovering
        assert(isequal(point.mask, oldPNNN.point.mask), 'Recovered masks differ.');
        checkCosts(point.row, oldPNNN.point.row);
        checkPrediction(point.fullSignalPrediction, oldPNNN.point.fullSignalPrediction, 1e-6);
    end
    models{3} = point.model;
    metrics = [linear.complexTable; linear.pniqTable; point.row];
    assert(all(metrics.ActualRealParameters == P), 'Active budget differs from P.');
    bundle = struct('version', 1, 'P', P, 'models', {models}, ...
        'config', cfg, 'signature', signature, 'split', split, ...
        'sourceFile', xyFile, 'sampleRateHz', source.fs, ...
        'trainingMeans', trainingMeans, 'metrics', metrics, ...
        'memoryConvention', "periodic whole record", ...
        'mapping', "desired u -> DC-removed ILC predistortion", ...
        'matlabVersion', version, 'rfValidationPerformed', false);
    inferred = applyDPDModels(bundle, u);
    checkPrediction(inferred(1).yvalmod, linear.predictions.complexFull, 1e-10);
    checkPrediction(inferred(2).yvalmod, linear.predictions.pniqFull, 1e-10);
    checkPrediction(inferred(3).yvalmod, point.fullSignalPrediction, 1e-6);
    save(modelFile, 'bundle', '-v7.3');
    reloaded = load(modelFile, 'bundle');
    roundTrip = applyDPDModels(reloaded.bundle, u);
    for k = 1:3
        checkPrediction(roundTrip(k).yvalmod, inferred(k).yvalmod, 1e-6);
    end
end

%% 5. Export fresh predictions, never overwrite old experimental files
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
outputFile = fullfile(outDir, ['experiment20260429T134032_xy_execution_' stamp '.mat']);
dpd = applyDPDModels(bundle, u, outputFile);
fprintf('\nOffline NMSE against DC-removed ILC target (NOT PA output NMSE):\n');
for k = 1:3
    fprintf('  %-24s %8.3f dB\n', dpd(k).modeltype, nmseComplexDb(y, dpd(k).yvalmod));
end
fprintf('Models: %s\nSignals: %s\nNo RF transmission performed.\n', modelFile, outputFile);
end

function checkPrediction(actual, expected, tolerance)
relativeError = norm(actual(:) - expected(:)) / max(norm(expected(:)), eps);
assert(relativeError <= tolerance, ...
    'Recovered/inferred prediction differs (relative error %.3g > %.3g).', ...
    relativeError, tolerance);
end

function checkCosts(actual, expected)
fields = {'ActualRealParameters','SelectedLambda','FLOPsPerSample','ActiveWeights','ActiveBiases'};
assert(isequaln(actual(:, fields), expected(:, fields)), ...
    'Recovered model changed the scientific budget, lambda or FLOPs.');
end
