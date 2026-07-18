function results = run_selected_comparison(selectedParameters, sweep)
% run_selected_comparison - Present one signed parameter-sweep point.
% The three main predictions come from checkpoints; six fixed-ridge
% predictions are refitted on the same stored supports before plotting.

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
for folder = ["metrics","pn_gmp_comparison","pnnn","splits","sweep"]
    addpath(fullfile(projectRoot, 'toolbox', folder));
end

if nargin < 2 || isempty(sweep)
    sweep = run_parameter_sweep(20:10:500);
end

%% Select the three comparable rows and load their signed checkpoints
identity = sweep.sweepIdentity;
signature = identity.experimentSignature;
fixedLambdas = identity.fixedRidgeLambdas(:);
if ~ismember(selectedParameters, identity.parameterGrid)
    error('run_selected_comparison:UnsignedTarget', ...
        'The requested target is not part of the signed grid.');
end
modelNames = ["Complex GMP DOMP sweep"; ...
    "Independent PN-IQ PN-DOMP sweep"; "Sparse PNNN N12"];
selectedRows = sweep.results([],:);
for model = modelNames.'
    row = string(sweep.results.Model) == model & ...
        sweep.results.TargetRealParameters == selectedParameters;
    if nnz(row) ~= 1
        error('run_selected_comparison:MissingSweepPoint', ...
            'Target %d must occur once in each model family.', ...
            selectedParameters);
    end
    selectedRows = [selectedRows; sweep.results(row,:)]; %#ok<AGROW>
end
if any(selectedRows.ActualRealParameters ~= selectedParameters)
    error('run_selected_comparison:ParameterMismatch', ...
        'Every selected model must have the requested real-parameter count.');
end

directory = char(string(sweep.resultDirectory));
stored = load(fullfile(directory, 'complexity_sweep.mat'), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The sweep summary belongs to another experiment.');
end
summaryPayload = artifact.payload;

stored = load(fullfile(directory, 'linear_sweep.mat'), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The linear checkpoint belongs to another experiment.');
end
linear = artifact.payload;

stored = load(fullfile(directory, 'sweep_dense_source.mat'), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The dense PNNN checkpoint belongs to another experiment.');
end
denseSource = artifact.payload;

stored = load(fullfile(directory, ...
    sprintf('pnnn_target_%04d.mat', selectedParameters)), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The sparse PNNN checkpoint belongs to another experiment.');
end
pnnn = artifact.payload;

if ~isequaln(summaryPayload.results, sweep.results) || ...
        ~isequaln(summaryPayload.evaluationProtocol, ...
        sweep.evaluationProtocol)
    error('run_selected_comparison:SummaryMismatch', ...
        'The in-memory sweep differs from its signed summary.');
end

complexIndex = find(linear.complexTable.TargetRealParameters == ...
        selectedParameters);
pnIndex = find(linear.pnTable.TargetRealParameters == selectedParameters);
if numel(complexIndex) ~= 1 || numel(pnIndex) ~= 1 || ...
        pnnn.row.TargetRealParameters ~= selectedParameters || ...
        pnnn.denseSourceDigest ~= denseSource.digest
    error('run_selected_comparison:ArtifactMismatch', ...
        'The selected linear or PNNN checkpoint is incompatible.');
end
if ~isequaln(linear.complexTable(complexIndex, :), selectedRows(1, :)) || ...
        ~isequaln(linear.pnTable(pnIndex, :), selectedRows(2, :)) || ...
        ~isequaln(pnnn.row, selectedRows(3, :))
    error('run_selected_comparison:ArtifactRowMismatch', ...
        'Artifact rows must equal the corresponding signed summary rows.');
end
complexPrediction = linear.predictions.complexFull(:, complexIndex);
pnPrediction = linear.predictions.pnFull(:, pnIndex);
pnnnPrediction = pnnn.fullSignalPrediction(:);

%% Recover the signed measurement and aligned full-signal target
cfg = getFairDOMPComparisonConfig(projectRoot);
measurement = load(cfg.measurementFile);
[x, y] = selectXYByMapping(measurement.x, measurement.y, cfg.mappingMode);
x = x(:);
y = y(:);
if ~isequaln(buildExperimentSignature(x, y, cfg), signature)
    error('run_selected_comparison:SignatureMismatch', ...
        'The configured measurement differs from the signed sweep.');
end
if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end
split = buildCommonComparisonSplit(x, y, cfg);
fullSignalIndices = double(split.fullSignalIndices(:));
targetFullSignal = y(fullSignalIndices);
sampleRateHz = resolveMeasurementSampleRate(measurement);

%% Refit the six fixed-ridge references on the selected supports
stored = load(fullfile(directory, ...
    'fixed_lambda_linear_sweep.mat'), 'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The fixed-ridge checkpoint belongs to another experiment.');
end
fixedSweep = artifact.payload;
selectedLinear = struct( ...
    'complexTable', linear.complexTable(complexIndex, :), ...
    'pnTable', linear.pnTable(pnIndex, :), ...
    'paths', linear.paths, ...
    'pnPathMap', linear.pnPathMap);
selectedConfig = cfg;
selectedConfig.sweep.parameterGrid = selectedParameters;
evaluated = run_fixed_ridge_sweep( ...
    x, y, split, selectedConfig, selectedLinear, true);
fixedTable = orderFixedRows( ...
    fixedSweep, selectedParameters, fixedLambdas);
evaluatedRows = orderFixedRows( ...
    evaluated.table, selectedParameters, fixedLambdas);
if any(evaluatedRows.ActualRealParameters ~= ...
        fixedTable.ActualRealParameters) || ...
        any(evaluatedRows.FLOPsPerSample ~= fixedTable.FLOPsPerSample) || ...
        any(abs(evaluatedRows.IdentificationNMSEdB - ...
        fixedTable.IdentificationNMSEdB) > 1e-9) || ...
        any(abs(evaluatedRows.FullSignalNMSEdB - ...
        fixedTable.FullSignalNMSEdB) > 1e-9)
    error('run_selected_comparison:FixedRidgeMismatch', ...
        'Selected ridge refits differ from the saved sweep.');
end
expectedFixedParameters = repelem( ...
    selectedRows.ActualRealParameters(1:2), numel(fixedLambdas));
expectedFixedFLOPs = repelem( ...
    selectedRows.FLOPsPerSample(1:2), numel(fixedLambdas));
if any(fixedTable.ActualRealParameters ~= expectedFixedParameters) || ...
        any(fixedTable.FLOPsPerSample ~= expectedFixedFLOPs)
    error('run_selected_comparison:FixedRidgeComplexityMismatch', ...
        'Fixed Ridge must retain the selected supports, parameters, and FLOPs.');
end
fixedPredictionMatrix = [evaluated.predictions.complexFull, ...
    evaluated.predictions.pnFull];

%% Build the public three-row comparison
comparisonTable = selectedRows;
comparisonTable.Model = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"; ...
    "Sparse PNNN N12"];
supportCount = selectedParameters/2;
pnFeatures = linear.paths.pn(1:supportCount);
results = struct('selectedParameters', selectedParameters, ...
    'comparisonTable', comparisonTable, ...
    'fixedLambdaComparisonTable', fixedTable, ...
    'linearSupports', struct( ...
        'complex', linear.paths.complex(1:supportCount), ...
        'pnFeatures', pnFeatures, ...
        'pnComplex', unique(linear.pnPathMap.SourceRegressorIndex( ...
            1:supportCount), 'stable')), ...
    'sampleRateHz', sampleRateHz, ...
    'evaluationProtocol', sweep.evaluationProtocol);

%% Compute one shared Welch spectrum and write four direct figures
mainPredictions = [complexPrediction, pnPrediction, pnnnPrediction];
spectrum = computeSelectedPointSpectra(targetFullSignal, mainPredictions, ...
    sampleRateHz, fixedPredictionMatrix);
selectedDirectory = fullfile(directory, ...
    sprintf('selected_point_%04d', selectedParameters));
if ~isfolder(selectedDirectory)
    mkdir(selectedDirectory);
end
outputFile = fullfile(selectedDirectory, 'selected_output_spectrum.png');
errorFile = fullfile(selectedDirectory, 'selected_error_spectrum.png');
ridgeOutputFile = fullfile(selectedDirectory, ...
    'selected_ridge_output_spectrum.png');
ridgeErrorFile = fullfile(selectedDirectory, ...
    'selected_ridge_error_spectrum.png');
plotCurves(spectrum.frequencyMHz, spectrum.outputPSDdB, ...
    ["Target full signal", "Complex GMP-DOMP", ...
    "PN-IQ PN-DOMP", "Sparse PNNN N12"], ...
    'Normalized PSD (dB)', outputFile);
plotCurves(spectrum.frequencyMHz, spectrum.errorPSDdB, ...
    ["Complex GMP-DOMP error", "PN-IQ PN-DOMP error", ...
    "Sparse PNNN N12 error"], ...
    'Error PSD relative to target peak (dB)', errorFile);
plotRidgePanels(spectrum.frequencyMHz, ...
    {[spectrum.outputPSDdB(:, [1 2]), spectrum.fixedOutputPSDdB(:, 1:3)], ...
    [spectrum.outputPSDdB(:, [1 3]), spectrum.fixedOutputPSDdB(:, 4:6)]}, ...
    true, 'Normalized PSD (dB)', ridgeOutputFile);
plotRidgePanels(spectrum.frequencyMHz, ...
    {[spectrum.errorPSDdB(:, 1), spectrum.fixedErrorPSDdB(:, 1:3)], ...
    [spectrum.errorPSDdB(:, 2), spectrum.fixedErrorPSDdB(:, 4:6)]}, ...
    false, 'Error PSD relative to target peak (dB)', ridgeErrorFile);

results.outputSpectrumFigure = string(outputFile);
results.errorSpectrumFigure = string(errorFile);
results.ridgeOutputSpectrumFigure = string(ridgeOutputFile);
results.ridgeErrorSpectrumFigure = string(ridgeErrorFile);
results.spectrumConfig = spectrum.config;
fprintf('\n=== Main comparison ===\n');
disp(comparisonTable);
fprintf('\n=== Fixed Ridge comparison ===\n');
disp(fixedTable);
fprintf('Reused linear, fixed-ridge, dense N12, and sparse point checkpoints: YES\n');
fprintf('Selected-point figures: %s\n', selectedDirectory);
end

function rows = orderFixedRows(value, target, fixedLambdas)
models = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"];
rows = value([],:);
for model = models.'
    for lambda = fixedLambdas.'
        match = string(value.Model) == model & ...
            value.TargetRealParameters == target & value.FixedLambda == lambda;
        if nnz(match) ~= 1
            error('run_selected_comparison:MissingFixedRow', ...
                'Each family and fixed lambda must occur exactly once.');
        end
        rows = [rows; value(match,:)]; %#ok<AGROW>
    end
end
end

function plotCurves(frequencyMHz, values, labels, yLabel, filename)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 940 600]);
cleanup = onCleanup(@() close(figureHandle));
lines = plot(frequencyMHz, values, 'LineWidth', 1.15);
for index = 1:numel(lines)
    lines(index).DisplayName = labels(index);
end
grid on;
xlabel('Frequency (MHz)');
ylabel(yLabel);
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', min(numel(lines), 4), 'Interpreter', 'none');
exportgraphics(figureHandle, filename, 'Resolution', 160, ...
    'BackgroundColor', 'white');
clear cleanup;
end

function plotRidgePanels(frequencyMHz, panels, includeTarget, yLabel, filename)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1120 780]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
familyNames = ["Complex GMP-DOMP", "PN-IQ PN-DOMP"];
familyColors = lines(2);
lineStyles = ["-", "--", "-.", ":"];
for panelIndex = 1:2
    axesHandle = nexttile(layout);
    hold(axesHandle, 'on');
    values = panels{panelIndex};
    if includeTarget
        plot(axesHandle, frequencyMHz, values(:, 1), ...
            'Color', [0.15 0.15 0.15], 'LineWidth', 1.2, ...
            'DisplayName', 'Target full signal');
        values = values(:, 2:end);
    end
    labels = [familyNames(panelIndex) + " principal", ...
        familyNames(panelIndex) + ", lambda=1e-3", ...
        familyNames(panelIndex) + ", lambda=1e-4", ...
        familyNames(panelIndex) + ", lambda=1e-5"];
    for lineIndex = 1:4
        plot(axesHandle, frequencyMHz, values(:, lineIndex), ...
            'Color', familyColors(panelIndex, :), ...
            'LineStyle', lineStyles(lineIndex), ...
            'LineWidth', 1 + (lineIndex == 1), ...
            'DisplayName', labels(lineIndex));
    end
    grid(axesHandle, 'on');
    xlabel(axesHandle, 'Frequency (MHz)');
    ylabel(axesHandle, yLabel);
    title(axesHandle, familyNames(panelIndex), 'Interpreter', 'none');
    legend(axesHandle, 'Location', 'eastoutside', 'Interpreter', 'none');
end
exportgraphics(figureHandle, filename, 'Resolution', 160, ...
    'BackgroundColor', 'white');
clear cleanup;
end
