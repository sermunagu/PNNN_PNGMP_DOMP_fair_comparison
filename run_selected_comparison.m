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
signature = sweep.experimentSignature;
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
if ~isequaln(artifact.experimentSignature, signature) || ...
        ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The sweep summary belongs to another experiment.');
end
summary = artifact.payload;

stored = load(fullfile(directory, 'linear_sweep.mat'), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.experimentSignature, signature) || ...
        ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The linear checkpoint belongs to another experiment.');
end
linear = artifact.payload;

stored = load(fullfile(directory, 'sweep_dense_source.mat'), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.experimentSignature, signature) || ...
        ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The dense PNNN checkpoint belongs to another experiment.');
end
denseSource = artifact.payload;

stored = load(fullfile(directory, ...
    sprintf('pnnn_target_%04d.mat', selectedParameters)), ...
    'checkpointArtifact');
artifact = stored.checkpointArtifact;
if ~isequaln(artifact.experimentSignature, signature) || ...
        ~isequaln(artifact.sweepIdentity, identity)
    error('run_selected_comparison:ExperimentMismatch', ...
        'The sparse PNNN checkpoint belongs to another experiment.');
end
pnnn = artifact.payload;

if ~isequaln(summary.results, sweep.results)
    error('run_selected_comparison:SummaryMismatch', ...
        'The in-memory sweep differs from its signed summary.');
end

complexIndex = find(linear.complexTable.TargetRealParameters == ...
    selectedParameters);
pnIndex = find(linear.pnTable.TargetRealParameters == selectedParameters);
if numel(complexIndex) ~= 1 || numel(pnIndex) ~= 1 || ...
        pnnn.target ~= selectedParameters || ...
        ~isequaln(pnnn.denseSourceSignature, denseSource.signature)
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

%% Recover the aligned full-signal target and fitting context
% Fixtures can provide this context directly; real runs verify the signature.
contextFields = {'targetFullSignal','fullSignalIndices','sampleRateHz'};
if all(isfield(sweep, contextFields))
    targetFullSignal = sweep.targetFullSignal(:);
    fullSignalIndices = double(sweep.fullSignalIndices(:));
    sampleRateHz = double(sweep.sampleRateHz);
    sampleRateSource = "Provided signed-sweep context";
    if isfield(sweep, 'sampleRateSource')
        sampleRateSource = string(sweep.sampleRateSource);
    end
    fittingContext = [];
else
    cfg = getFairDOMPComparisonConfig(projectRoot);
    measurement = load(cfg.measurementFile, 'x', 'y', 'fs', 'info_signal');
    [x, y] = selectXYByMapping( ...
        measurement.x, measurement.y, cfg.mappingMode);
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
    % TODO: review the scientific split protocol after this equivalence rewrite.
    fullSignalIndices = double(split.fullSignalIndices(:));
    targetFullSignal = y(fullSignalIndices);
    [sampleRateHz, sampleRateSource] = measurementSampleRate(measurement);
    fittingContext = struct('x', x, 'y', y, 'split', split, 'cfg', cfg);
end
%% Refit the six fixed-ridge references on the selected supports
providedFixed = {'fixedLambdaComparisonTable', ...
    'fixedLambdaFullSignalPredictions'};
if all(isfield(sweep, providedFixed))
    fixedTable = sweep.fixedLambdaComparisonTable;
    fixedPredictions = sweep.fixedLambdaFullSignalPredictions;
    ridgeTime = 0;
else
    stored = load(fullfile(directory, ...
        'fixed_lambda_linear_sweep.mat'), 'checkpointArtifact');
    artifact = stored.checkpointArtifact;
    if ~isequaln(artifact.experimentSignature, signature) || ...
            ~isequaln(artifact.sweepIdentity, identity)
        error('run_selected_comparison:ExperimentMismatch', ...
            'The fixed-ridge checkpoint belongs to another experiment.');
    end
    fixedSweep = artifact.payload;
    selectedLinear = linear;
    selectedLinear.complexTable = linear.complexTable(complexIndex, :);
    selectedLinear.pnTable = linear.pnTable(pnIndex, :);
    selectedLinear.supports.complex = linear.supports.complex(complexIndex);
    selectedLinear.supports.pnFeatures = ...
        linear.supports.pnFeatures(pnIndex);
    selectedLinear.supports.pnComplex = linear.supports.pnComplex(pnIndex);
    selectedConfig = fittingContext.cfg;
    selectedConfig.sweep.parameterGrid = selectedParameters;
    timer = tic;
    evaluated = run_fixed_ridge_sweep(fittingContext.x, ...
        fittingContext.y, fittingContext.split, selectedConfig, ...
        selectedLinear, true);
    ridgeTime = toc(timer);
    fixedTable = orderFixedRows(fixedSweep.table, selectedParameters);
    evaluatedRows = orderFixedRows(evaluated.table, selectedParameters);
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
    fixedPredictions = struct( ...
        'complexGMP', predictionFamily(evaluated.predictions.complexFull), ...
        'pnIQ', predictionFamily(evaluated.predictions.pnFull));
end
fixedTable = orderFixedRows(fixedTable, selectedParameters);
expectedFixedParameters = repelem(selectedRows.ActualRealParameters(1:2), 3);
expectedFixedFLOPs = repelem(selectedRows.FLOPsPerSample(1:2), 3);
if any(fixedTable.ActualRealParameters ~= expectedFixedParameters) || ...
        any(fixedTable.FLOPsPerSample ~= expectedFixedFLOPs)
    error('run_selected_comparison:FixedRidgeComplexityMismatch', ...
        'Fixed Ridge must retain the selected supports, parameters, and FLOPs.');
end
fixedPredictionMatrix = [fixedPredictions.complexGMP.lambda1e3(:), ...
    fixedPredictions.complexGMP.lambda1e4(:), ...
    fixedPredictions.complexGMP.lambda1e5(:), ...
    fixedPredictions.pnIQ.lambda1e3(:), ...
    fixedPredictions.pnIQ.lambda1e4(:), ...
    fixedPredictions.pnIQ.lambda1e5(:)];

%% Build the public three-row comparison
comparisonTable = selectedRows(:, {'Model','ActualRealParameters', ...
    'FullSignalNMSEdB','FLOPsPerSample'});
comparisonTable.Model = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"; ...
    "Sparse PNNN N12"];
results = struct('selectedParameters', selectedParameters, ...
    'comparisonTable', comparisonTable, ...
    'fixedLambdaComparisonTable', fixedTable, ...
    'resultDirectory', string(directory), ...
    'sweepResultDirectory', string(directory), ...
    'reusedLinearSweep', true, 'reusedPNNNPoint', true, ...
    'reusedDenseSource', true, 'selectedSweepRows', selectedRows, ...
    'linearSupports', struct( ...
        'complex', linear.supports.complex{complexIndex}, ...
        'pnFeatures', linear.supports.pnFeatures{pnIndex}, ...
        'pnComplex', linear.supports.pnComplex{pnIndex}), ...
    'targetFullSignal', targetFullSignal, ...
    'fullSignalIndices', fullSignalIndices, ...
    'sampleRateHz', sampleRateHz, ...
    'sampleRateSource', sampleRateSource, ...
    'ridgePredictionTimeSeconds', ridgeTime, ...
    'fullSignalPredictions', struct('complexGMP', complexPrediction, ...
        'pnIQ', pnPrediction, 'sparsePNNNN12', pnnnPrediction), ...
    'fixedLambdaFullSignalPredictions', fixedPredictions);

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

results.selectedPointDirectory = string(selectedDirectory);
results.outputSpectrumFigure = string(outputFile);
results.errorSpectrumFigure = string(errorFile);
results.ridgeOutputSpectrumFigure = string(ridgeOutputFile);
results.ridgeErrorSpectrumFigure = string(ridgeErrorFile);
results.spectrumConfig = spectrum.config;
results.spectrumConfig.sampleRateSource = sampleRateSource;
disp(comparisonTable);
fprintf('Reused linear, fixed-ridge, dense N12, and sparse point checkpoints: YES\n');
fprintf('Selected fixed-ridge predictions completed in %.2f s.\n', ridgeTime);
fprintf('Selected-point figures: %s\n', selectedDirectory);
end

function rows = orderFixedRows(value, target)
models = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"];
lambdas = [1e-3; 1e-4; 1e-5];
rows = value([],:);
for model = models.'
    for lambda = lambdas.'
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

function predictions = predictionFamily(values)
predictions = struct('lambda1e3', values(:, 1), ...
    'lambda1e4', values(:, 2), 'lambda1e5', values(:, 3));
end

function [sampleRateHz, source] = measurementSampleRate(measurement)
hasFs = isfield(measurement, 'fs');
hasInfoFs = isfield(measurement, 'info_signal') && ...
    isstruct(measurement.info_signal) && ...
    isfield(measurement.info_signal, 'fsovs');
if hasFs
    sampleRateHz = double(measurement.fs);
    source = "measurement.fs";
    if hasInfoFs
        infoSampleRate = double(measurement.info_signal.fsovs);
        tolerance = 1e-9*max([1, abs(sampleRateHz), abs(infoSampleRate)]);
        if abs(infoSampleRate - sampleRateHz) > tolerance
            error('run_selected_comparison:SampleRateMismatch', ...
                'Measurement fs and info_signal.fsovs disagree.');
        end
        source = "measurement.fs, confirmed by info_signal.fsovs";
    end
else
    sampleRateHz = double(measurement.info_signal.fsovs);
    source = "measurement.info_signal.fsovs";
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
