function run_simplification_regression_test(referenceFile, recordReference)
% Record once before refactoring; compare numerical outputs after each block.
if nargin < 2, recordReference = false; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'), genpath(fullfile(root, 'toolbox')));
cfg = getFairDOMPComparisonConfig(root);
cfg.sweep.parameterGrid = [20 30];
cfg.identificationFraction = 0.75;
cfg.training.miniBatchSize = 32;
cfg.training.historicalMaxEpochs = 2;
cfg.training.historicalLearnRateDropPeriod = 1;
cfg.training.historicalValidationPatience = 2;
cfg.training.verbose = false;
cfg.pruning.historicalFineTuneEpochs = 2;
cfg.pnnn.sparseBaseHiddenNeurons = 3;
rng(715, 'twister');
x = complex(randn(160, 1), randn(160, 1))/sqrt(2);
y = 0.8*x + 0.12*circshift(x, 1).*abs(circshift(x, 2)).^2;
split = buildCommonComparisonSplit(x, y, cfg);
manager = GMP_createRegressorManager(x, y, cfg.gmp);
actual.population = compactGMPRegressors(manager, 1:numel(manager.regPopulation));
assert(numel(actual.population) == 673);
[actual.linear, actual.models] = run_linear_sweep(x, y, split, cfg);
actual.fixed = run_fixed_ridge_sweep(x, y, split, cfg, actual.linear, true);
[source, features, targets, rotation] = prepare_pnnn_dense_source( ...
    x, y, split, cfg, cfg.reducedRealParameterTarget);
source.runtimeConfig.training.verbose = false;
actual.denseWeights = source.denseFit.network.Learnables;
actual.normalization = source.denseFit.normalization;
actual.epochs = [source.bestDenseEpoch source.fineTuneEpochs];
for k = 1:2
    point = fit_sparse_pnnn_target(source, cfg.sweep.parameterGrid(k), ...
        features, targets, rotation, y, split, cfg);
    actual.points{k} = struct('row', point.row, 'mask', {point.mask}, ...
        'prediction', point.fullSignalPrediction, ...
        'weights', point.model.network.Learnables);
end
assert(isequaln(source.denseFit.network.Learnables, actual.denseWeights));
if recordReference
    assert(~isfile(referenceFile), 'Do not overwrite the pre-refactor reference.');
    save(referenceFile, 'actual', '-v7.3');
    fprintf('SIMPLIFICATION REFERENCE SAVED: %s\n', referenceFile);
else
    before = load(referenceFile, 'actual');
    names = fieldnames(actual);
    for k = 1:numel(names)
        assert(isequaln(actual.(names{k}), before.actual.(names{k})), ...
            'Pre/post difference in %s', names{k});
    end
    fprintf('SIMPLIFICATION REGRESSION: PASS (exact equality, all fields)\n');
end
end
