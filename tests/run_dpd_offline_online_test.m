function run_dpd_offline_online_test(folder)
% Two tiny fresh fits; verify model budgets, reload and inference without training.
% Pass a folder to keep the fixture for the fresh-session test.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'toolbox')));
if nargin < 1
    folder = tempname;
    cleanup = onCleanup(@() rmdir(folder, 's'));
end
mkdir(folder);
cfg = getFairDOMPComparisonConfig(root);
cfg.training.historicalMaxEpochs = 2;
cfg.training.historicalLearnRateDropPeriod = 1;
cfg.training.historicalValidationPatience = 2;
cfg.pruning.historicalFineTuneEpochs = 2;
cfg.dpd.resultsRoot = fullfile(folder, 'offline');
rng(932, 'twister');
x = 0.3*complex(randn(512,1), randn(512,1));
y = 0.8*x + 0.09*circshift(x,1).*abs(circshift(x,2)).^2;
inputFile = fullfile(folder, 'input.mat');
save(inputFile, 'x');
modelFiles = cell(1, 2);
profile clear; profile on;
for k = 1:2
    measurementFile = fullfile(folder, sprintf('capture%d.mat', k));
    save(measurementFile, 'x', 'y');
    [modelFiles{k}, bundle] = trainDPDOffline(measurementFile, 20, cfg);
    assert(all(bundle.metrics.ActualRealParameters == 20));
    assert(2*numel(bundle.models{1}.coefficients) == 20);
    assert(numel(bundle.models{2}.coefficientsI) + numel(bundle.models{2}.coefficientsQ) == 20);
    assert(sum(cellfun(@nnz, bundle.models{3}.masks)) == 20);
    x = 1.3*x;
    y = 0.6*x + 0.2j*circshift(x,2).*abs(x).^2;
end
profile off; info = profile('info');
functions = info.FunctionTable;
assert(sum([functions(strcmp({functions.FunctionName}, 'selectDOMPSupport')).NumCalls]) == 4);
assert(sum([functions(strcmp({functions.FunctionName}, 'trainnet')).NumCalls]) == 4);
assert(~strcmp(modelFiles{1}, modelFiles{2}));
modelA = modelFiles{1};
direct = load(fullfile(fileparts(modelA), 'predictions.mat'));
profile clear; profile on;
[outputFile, expected] = runDPDOnline(modelA, inputFile);
secondOutput = runDPDOnline(modelA, inputFile);
profile off; info = profile('info');
for name = ["selectDOMPSupport", "prepare_pnnn_dense_source", "pruneAndFineTunePNNN", "trainnet"]
    assert(~any(strcmp({info.FunctionTable.FunctionName}, name)));
end
assert(~strcmp(outputFile, secondOutput));
assert(isequal([expected.yvalmod], direct.predictions));
saved = load(outputFile);
assert(isequal(saved.dpd, expected));
assert(isequal(size([expected.yvalmod]), [512 3]));
assert(~isreal([expected.yvalmod]) && all(isfinite([expected.yvalmod]), 'all'));
save(fullfile(folder, 'fresh_session_fixture.mat'), 'modelA', 'inputFile', 'expected');
fprintf('DPD OFFLINE/ONLINE: PASS.\n');
end
