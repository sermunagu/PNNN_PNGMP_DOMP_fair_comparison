function sweep = runDPDSweep(measurementFile, parameterGrid, cfg)
% Sweep u -> x_ILC, or resume explicitly by passing an existing run folder.
% Optional cfg supports small tests; scientific fits are shared with OFFLINE.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'config'), genpath(fullfile(root, 'toolbox')));
if isfolder(measurementFile)
    directory = char(measurementFile);
    saved = load(fullfile(directory, 'experiment.mat'));
    cfg = saved.cfg;
    x = saved.x;
    y = saved.y;
    split = saved.split;
    provenance = saved.provenance;
else
    if nargin < 3, cfg = getComparisonConfig(root); end
    if nargin < 2, parameterGrid = 20:10:500; end
    cfg.measurementFile = measurementFile;
    cfg.sweep.parameterGrid = unique(double(parameterGrid(:).'));
    if ~isfield(cfg, 'dpd'), cfg.dpd = struct(); end
    if ~isfield(cfg.dpd, 'sweepResultsRoot')
        cfg.dpd.sweepResultsRoot = fullfile(root, 'results', 'dpd_sweep');
    end
    if ~isfield(cfg.dpd, 'exportFigures'), cfg.dpd.exportFigures = true; end
    capture = load(measurementFile);
    x = double(capture.x(:));
    y = double(capture.y(:));
    provenance = struct('sourceFile', measurementFile, ...
        'inputMean', mean(x), 'outputMean', mean(y), ...
        'signalDomain', "x = u; y = x_ILC; offline ILC imitation", ...
        'rfValidationPerformed', false);
    if isfield(capture, 'fs'), provenance.fs = capture.fs; end
    if cfg.pnnn.removeDC
        x = x - mean(x);
        y = y - mean(y);
    end
    split = buildCommonComparisonSplit(x, y, cfg);
    % Historical target-informed split is unchanged; this is not an external test.
    split.evaluationIndices = setdiff(split.fullSignalIndices, split.identificationIndices);
    if ~isfolder(cfg.dpd.sweepResultsRoot), mkdir(cfg.dpd.sweepResultsRoot); end
    directory = tempname(cfg.dpd.sweepResultsRoot);
    mkdir(directory);
    save(fullfile(directory, 'experiment.mat'), 'cfg', 'x', 'y', 'split', 'provenance', '-v7.3');
end
fprintf('DPD sweep directory (use this path to resume): %s\n', directory);
targets = cfg.sweep.parameterGrid;
evaluationRows = split.evaluationIndices;

%% One DOMP path per linear family, with all prefix fits and inferable models
filename = fullfile(directory, 'linear_sweep.mat');
if isfile(filename)
    saved = load(filename, 'linear', 'linearModels');
    linear = saved.linear;
    linearModels = saved.linearModels;
else
    [linear, linearModels] = fitLinearModels(x, y, split, cfg);
    save(filename, 'linear', 'linearModels', '-v7.3');
end
filename = fullfile(directory, 'fixed_lambda_linear_sweep.mat');
if isfile(filename)
    saved = load(filename, 'fixedResults');
    fixedResults = saved.fixedResults;
else
    fixed = fitFixedRidge(x, y, split, cfg, linear, true);
    fixedResults = fixed.table;
    predictions = [fixed.predictions.complexFull, fixed.predictions.pniqFull];
    for index = 1:height(fixedResults)
        fixedResults.EvaluationNMSEdB(index) = nmseComplexDb( ...
            y(evaluationRows), predictions(evaluationRows, index));
    end
    save(filename, 'fixedResults', '-v7.3');
    clear fixed predictions;
end

%% Shared definitive dense network; each P starts from this same source
pointFiles = fullfile(directory, compose("pnnn_target_%04d.mat", targets));
if any(~isfile(pointFiles))
    filename = fullfile(directory, 'sweep_dense_source.mat');
    if isfile(filename)
        saved = load(filename, 'denseSource');
        denseSource = saved.denseSource;
        [features, neuralTargets, rotation] = buildPhaseNormDataset( ...
            x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode, cfg.pnnn);
        features = features.';
        neuralTargets = neuralTargets.';
        rotation = rotation(:);
    else
        [denseSource, features, neuralTargets, rotation] = trainDensePNNN( ...
            x, y, split, cfg, cfg.reducedRealParameterTarget);
        save(filename, 'denseSource', '-v7.3');
    end
    denseSource.runtimeConfig.training.verbose = cfg.training.verbose;
end
results = table();
modelFiles = strings(numel(targets), 1);
for index = 1:numel(targets)
    P = targets(index);
    filename = pointFiles(index);
    if isfile(filename)
        saved = load(filename, 'point');
        point = saved.point;
    else
        point = fitSparsePNNN(denseSource, P, features, ...
            neuralTargets, rotation, y, split, cfg);
        save(filename, 'point', '-v7.3');
    end
    metrics = [linear.complexTable(index,:); linear.pniqTable(index,:); point.row];
    % Keep families separate: concatenating with single PNNN output casts doubles.
    predictions = {linear.predictions.complexFull(:,index), ...
        linear.predictions.pniqFull(:,index), point.fullSignalPrediction};
    for family = 1:3
        metrics.EvaluationNMSEdB(family) = nmseComplexDb( ...
            y(evaluationRows), predictions{family}(evaluationRows));
    end
    results = [results; metrics]; %#ok<AGROW>
    modelFiles(index) = fullfile(directory, sprintf('models_P%04d.mat', P));
    if ~isfile(modelFiles(index))
        models = compactLinearModels(linearModels, index, P);
        models{3} = point.model;
        pointConfig = cfg;
        pointConfig.sweep.parameterGrid = P;
        bundle = struct('P', P, 'models', {models}, 'config', pointConfig, ...
            'sourceFile', provenance.sourceFile, 'split', split, ...
            'metrics', metrics, 'provenance', provenance);
        if isfield(provenance, 'fs'), bundle.fs = provenance.fs; end
        save(modelFiles(index), 'bundle', '-v7.3');
    end
end

%% Selection uses evaluation only and requires complete window coverage
selection = [];
selectionInput = results;
selectionInput.FullSignalNMSEdB = results.EvaluationNMSEdB;
try
    selection = selectOperatingPoint(selectionInput, cfg.selection);
    selectionReason = selection.summarySentence;
    writetable(selection.diagnosticsTable, fullfile(directory, 'operating_point_selection.csv'));
    writetable(selection.sensitivityTable, fullfile(directory, 'operating_point_selection_sensitivity.csv'));
catch exception
    if ~ismember(string(exception.identifier), [ ...
            "selectOperatingPoint:NoJointlyStabilizedPoint", ...
            "selectOperatingPoint:InsufficientWindowCoverage", ...
            "selectOperatingPoint:NonmonotonicFLOPs"])
        rethrow(exception);
    end
    selectionReason = "No operating point selected: " + string(exception.message);
    % A resumed sparse pilot must not leave an obsolete selected-P CSV behind.
    for name = ["operating_point_selection.csv", "operating_point_selection_sensitivity.csv"]
        selectionFile = fullfile(directory, name);
        if isfile(selectionFile)
            previous = readtable(selectionFile);
            writetable(previous([],:), selectionFile);
        end
    end
end
fprintf('%s\n', selectionReason);
evaluationProtocol = struct('selectionMetric', "EvaluationNMSEdB", ...
    'evaluationIsIndependentTest', false, ...
    'description', "Complement of historical identification rows; not used for fitting. " + ...
        "Target-informed segment selection, full-record DC and periodic memory remain. " + ...
        "Selecting P on these curves makes evaluation a development set, not independent confirmation.");
coefficientMetadata = struct('definition', cfg.sweep.coefficientRangeDefinition, ...
    'description', "Linear: unit-peak output and per-column peak-normalized regressors; " + ...
        "complex modulus for GMP, joint absolute I/Q maximum for PN-IQ. " + ...
        "PNNN: all active weights/biases in z-score parameterization; not identical physical gains.");
writetable(results, fullfile(directory, 'dpd_complexity_sweep.csv'));
writetable(fixedResults, fullfile(directory, 'fixed_lambda_linear_sweep.csv'));
% The investigator can fill this CSV later; resuming must preserve RF annotations.
rfFile = fullfile(directory, 'rf_measurements.csv');
if ~isfile(rfFile)
    rf = table(results.Model, results.TargetRealParameters, repelem(modelFiles, 3, 1), ...
        repmat(string(provenance.sourceFile), height(results), 1), strings(height(results), 1), ...
        'VariableNames', {'Model','P','ModelFile','CaptureFile','RFMeasurementFile'});
    writetable(rf, rfFile);
end
sweep = struct('results', results, 'fixedResults', fixedResults, 'config', cfg, ...
    'resultDirectory', string(directory), 'modelFiles', modelFiles, ...
    'provenance', provenance, 'evaluationProtocol', evaluationProtocol, ...
    'coefficientMetadata', coefficientMetadata, 'selection', selection, ...
    'selectionReason', selectionReason, 'rfManifestFile', string(rfFile));
filename = fullfile(directory, 'dpd_complexity_sweep.mat');
save(filename, 'sweep', '-v7.3');
if cfg.dpd.exportFigures
    sweep.figureFiles = plotDPDSweepFigures(sweep);
    save(filename, 'sweep', '-v7.3');
end
end

function models = compactLinearModels(source, index, P)
% Store only the active prefix, not the max-budget dictionary with zero weights.
count = P/2;
gmp = source{1};
gmp.regPopulation = gmp.regPopulation(1:count);
gmp.support = gmp.support(1:count);
gmp.coefficients = gmp.coefficients(1:count, index);
gmp.comparisonCoefficients = gmp.comparisonCoefficients(1:count, index);
pn = source{2};
allSources = unique(pn.featureMap.SourceRegressorIndex, 'stable');
pn.featureMap = pn.featureMap(1:count,:);
activeSources = unique(pn.featureMap.SourceRegressorIndex, 'stable');
[~, positions] = ismember(activeSources, allSources);
pn.regPopulation = pn.regPopulation(positions);
pn.descriptors = pn.descriptors(positions);
pn.support = pn.support(1:count);
[~, pn.selectedColumns] = ismember(pn.featureMap.SourceRegressorIndex, activeSources);
pn.selectedColumns(pn.featureMap.IsQ) = pn.selectedColumns(pn.featureMap.IsQ) + numel(activeSources);
pn.coefficientsI = pn.coefficientsI(1:count, index);
pn.coefficientsQ = pn.coefficientsQ(1:count, index);
pn.comparisonI = pn.comparisonI(1:count, index);
pn.comparisonQ = pn.comparisonQ(1:count, index);
models = {gmp, pn};
end
