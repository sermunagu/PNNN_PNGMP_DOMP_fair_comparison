function result = runPASweep(parameterGrid, resultsRoot, cfg)
% Run a new PA sweep, or resume explicitly with runPASweep(runDirectory).
% The optional cfg is useful for small synthetic experiments.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'), genpath(fullfile(projectRoot, 'toolbox')));
if nargin >= 1 && (ischar(parameterGrid) || isstring(parameterGrid))
    resultDirectory = char(parameterGrid);
    saved = load(fullfile(resultDirectory, 'experiment.mat'));
    cfg = saved.cfg;
    x = saved.x;
    y = saved.y;
    split = saved.split;
else
    if nargin < 3, cfg = getComparisonConfig(projectRoot); end
    if nargin >= 2 && ~isempty(resultsRoot), cfg.sweep.resultsRoot = resultsRoot; end
    if nargin >= 1 && ~isempty(parameterGrid)
        cfg.sweep.parameterGrid = unique(double(parameterGrid(:).'));
    end
    measurement = load(cfg.measurementFile, 'x', 'y');
    [x, y] = mapModelSignals(measurement.x, measurement.y, cfg.mappingMode);
    if cfg.pnnn.removeDC
        x = x - mean(x);
        y = y - mean(y);
    end
    split = buildCommonComparisonSplit(x, y, cfg);
    if ~isfolder(cfg.sweep.resultsRoot), mkdir(cfg.sweep.resultsRoot); end
    resultDirectory = tempname(cfg.sweep.resultsRoot);
    mkdir(resultDirectory);
    % Store the data used by this run so resuming cannot select another capture.
    save(fullfile(resultDirectory, 'experiment.mat'), 'cfg', 'x', 'y', 'split', '-v7.3');
end
targets = cfg.sweep.parameterGrid;
fprintf('Sweep directory (use this path to resume): %s\n', resultDirectory);

%% Linear DOMP paths and fixed Ridge
linearFile = fullfile(resultDirectory, 'linear_sweep.mat');
if isfile(linearFile)
    saved = load(linearFile, 'linear');
    linear = saved.linear;
else
    linear = fitLinearModels(x, y, split, cfg);
    save(linearFile, 'linear', '-v7.3');
end
fixedFile = fullfile(resultDirectory, 'fixed_lambda_linear_sweep.mat');
if isfile(fixedFile)
    saved = load(fixedFile, 'fixedResults');
    fixedResults = saved.fixedResults;
else
    fixed = fitFixedRidge(x, y, split, cfg, linear);
    fixedResults = fixed.table;
    save(fixedFile, 'fixedResults', '-v7.3');
end

%% One shared dense source; each sparse target starts independently from it
pointFiles = fullfile(resultDirectory, compose("pnnn_target_%04d.mat", targets));
if any(~isfile(pointFiles))
    denseFile = fullfile(resultDirectory, 'sweep_dense_source.mat');
    if isfile(denseFile)
        saved = load(denseFile, 'denseSource');
        denseSource = saved.denseSource;
        [features, neuralTargets, phaseRotation] = buildPhaseNormDataset( ...
            x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode, cfg.pnnn);
        features = features.';
        neuralTargets = neuralTargets.';
        phaseRotation = phaseRotation(:);
    else
        [denseSource, features, neuralTargets, phaseRotation] = ...
            trainDensePNNN(x, y, split, cfg, cfg.reducedRealParameterTarget);
        save(denseFile, 'denseSource', '-v7.3');
    end
    denseSource.runtimeConfig.training.verbose = cfg.training.verbose;
end
pnnnRows = linear.complexTable([],:);
for index = 1:numel(targets)
    filename = pointFiles(index);
    if isfile(filename)
        saved = load(filename, 'point');
        point = saved.point;
    else
        point = fitSparsePNNN(denseSource, targets(index), ...
            features, neuralTargets, phaseRotation, y, split, cfg);
        save(filename, 'point', '-v7.3');
    end
    pnnnRows(index,:) = point.row;
end

%% Metrics and scientific selection
results = [linear.complexTable; linear.pniqTable; pnnnRows];
evaluationProtocol = struct( ...
    'internalValidationUsedForSelectionOnly', true, ...
    'internalValidationConsumer', "PNNN only", ...
    'linearIdentificationScope', string(cfg.sweep.linearIdentificationScope), ...
    'linearPrincipalLambda', double(cfg.sweep.linearPrincipalLambda), ...
    'linearLambdaSelection', string(cfg.sweep.linearLambdaSelection), ...
    'fixedRidgeSupportPolicy', string(cfg.sweep.fixedRidgeSupportPolicy), ...
    'completeSignalIncludesIdentification', all(ismember(split.identificationIndices, split.fullSignalIndices)), ...
    'completeSignalIsIndependentHoldout', ~any(ismember(split.identificationIndices, split.fullSignalIndices)));

selection = selectOperatingPoint(results, cfg.selection);
coefficientMetadata = struct( ...
    'definition', string(cfg.sweep.coefficientRangeDefinition), ...
    'description', "Maximum absolute coefficient/parameter. Complex GMP " + ...
        "uses the complex modulus of coefficients fitted on selected " + ...
        "per-column peak-normalized regressors against unit-peak output; " + ...
        "PN-IQ-GMP uses " + ...
        "the joint absolute maximum of its I/Q coefficients fitted on " + ...
        "selected per-column peak-normalized real features against " + ...
        "phase-normalized, " + ...
        "unit-peak output; " + ...
        "PNNN parameters remain in their z-score-normalized parameterization.", ...
    'warning', "The metric compares numerical dynamic range for " + ...
        "implementation and quantization. Linear families use unit-peak " + ...
        "output and per-column peak-normalized regressors whereas PNNN " + ...
        "uses its z-score parameterization; magnitudes must not be " + ...
        "interpreted as identical physical gains.");
warning('runPASweep:ParameterizationWarning', '%s', ...
    coefficientMetadata.warning);
fprintf('%s\n', selection.summarySentence);

writetable(results, fullfile(resultDirectory, 'complexity_sweep.csv'));
writetable(fixedResults, fullfile(resultDirectory, 'fixed_lambda_linear_sweep.csv'));
writetable(selection.diagnosticsTable, fullfile(resultDirectory, 'operating_point_selection.csv'));
writetable(selection.sensitivityTable, fullfile(resultDirectory, 'operating_point_selection_sensitivity.csv'));
writelines([selection.summarySentence; coefficientMetadata.warning], ...
    fullfile(resultDirectory, 'operating_point_selection_summary.txt'));
result = struct('results', results, 'config', cfg, ...
    'evaluationProtocol', evaluationProtocol, 'selection', selection, ...
    'coefficientMetadata', coefficientMetadata, 'resultDirectory', string(resultDirectory));
summaryFile = fullfile(resultDirectory, 'complexity_sweep.mat');
save(summaryFile, 'result', '-v7.3');

%% Export only after numerical results have been saved
addpath(cfg.paper.matlab2tikzSource);
exportOptions = struct();
nmseOptions = struct('metricVariable', 'FullSignalNMSEdB', ...
    'metricLabel', 'Full-signal NMSE (dB)', ...
    'names', cfg.names, 'includeFixed', true, ...
    'fixedLambdas', cfg.fixedRidgeLambdas, 'isNMSE', true, ...
    'annotateSelected', false, 'selection', selection, ...
    'exportOptions', exportOptions);
fprintf('[Output] Writing NMSE-vs-parameters figure...\n');
figureFiles.nmseParameters = plotSweepPaperFigure(results, fixedResults, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(resultDirectory, 'comparison_nmse_parameters_sweep'), ...
    nmseOptions);
fprintf('[Output] Writing NMSE-vs-FLOPs figure...\n');
nmseOptions.includeFixed = false;
figureFiles.nmseFLOPs = plotSweepPaperFigure(results, table(), ...
    'FLOPsPerSample', 'FLOPs per sample', ...
    fullfile(resultDirectory, 'comparison_nmse_flops_sweep'), ...
    nmseOptions);
fprintf('[Output] Writing coefficient-range figure...\n');
rangeOptions = struct('metricVariable', 'MaxAbsRealParameter', ...
    'metricLabel', 'Maximum absolute coefficient/parameter', ...
    'includeFixed', true, 'fixedLambdas', cfg.fixedRidgeLambdas, ...
    'names', cfg.names, 'useLogWhenPositive', true, ...
    'annotateSelected', false, ...
    'exportOptions', exportOptions);
figureFiles.maxAbsParameter = plotSweepPaperFigure(results, fixedResults, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(resultDirectory, 'comparison_max_abs_parameter_sweep'), ...
    rangeOptions);

result.figureFiles = figureFiles;
save(summaryFile, 'result', '-v7.3');
disp(results(:, {'Model','ActualRealParameters','SelectedLambda', ...
    'FullSignalNMSEdB','FLOPsPerSample','MaxAbsRealParameter'}));
end
