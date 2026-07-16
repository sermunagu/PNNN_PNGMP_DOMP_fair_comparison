function results = writeSelectedPointSpectra(results)
% writeSelectedPointSpectra - Plot output and error spectra for one point.
% Main and fixed-ridge predictions share one aligned full-signal target PSD.
% Four lightweight PNG figures are written inside the selected-point folder.

required = {'selectedParameters','resultDirectory','targetFullSignal', ...
    'fullSignalIndices','sampleRateHz','sampleRateSource', ...
    'fullSignalPredictions','fixedLambdaFullSignalPredictions', ...
    'fixedLambdaComparisonTable','comparisonTable'};
if ~isstruct(results) || ~all(isfield(results, required)) || ...
        ~istable(results.comparisonTable) || height(results.comparisonTable) ~= 3
    error('writeSelectedPointSpectra:InvalidResults', ...
        'A complete three-family selected-point result is required.');
end
predictionFields = {'complexGMP','pnIQ','sparsePNNNN12'};
if ~all(isfield(results.fullSignalPredictions, predictionFields))
    error('writeSelectedPointSpectra:MissingPredictions', ...
        'The three selected full-signal predictions are required.');
end
predictions = [results.fullSignalPredictions.complexGMP(:), ...
    results.fullSignalPredictions.pnIQ(:), ...
    results.fullSignalPredictions.sparsePNNNN12(:)];
fixedPredictions = fixedPredictionMatrix( ...
    results.fixedLambdaFullSignalPredictions);
spectrum = computeSelectedPointSpectra( ...
    results.targetFullSignal, predictions, results.sampleRateHz, ...
    fixedPredictions);

directory = fullfile(char(string(results.resultDirectory)), ...
    sprintf('selected_point_%04d', results.selectedParameters));
if ~isfolder(directory)
    mkdir(directory);
end
outputFile = fullfile(directory, 'selected_output_spectrum.png');
errorFile = fullfile(directory, 'selected_error_spectrum.png');
ridgeOutputFile = fullfile(directory, ...
    'selected_ridge_output_spectrum.png');
ridgeErrorFile = fullfile(directory, ...
    'selected_ridge_error_spectrum.png');
outputLabels = ["Target full signal", "Complex GMP-DOMP", ...
    "PN-IQ PN-DOMP", "Sparse PNNN N12"];
errorLabels = ["Complex GMP-DOMP error", "PN-IQ PN-DOMP error", ...
    "Sparse PNNN N12 error"];
outputCurveCount = writeSpectrumFigure(spectrum.frequencyMHz, ...
    spectrum.outputPSDdB, outputLabels, 'Normalized PSD (dB)', outputFile);
errorCurveCount = writeSpectrumFigure(spectrum.frequencyMHz, ...
    spectrum.errorPSDdB, errorLabels, ...
    'Error PSD relative to target peak (dB)', errorFile);
ridgeOutputCounts = writeRidgeSpectrumFigure(spectrum.frequencyMHz, ...
    {[spectrum.outputPSDdB(:, [1 2]), spectrum.fixedOutputPSDdB(:, 1:3)], ...
    [spectrum.outputPSDdB(:, [1 3]), spectrum.fixedOutputPSDdB(:, 4:6)]}, ...
    true, 'Normalized PSD (dB)', ridgeOutputFile);
ridgeErrorCounts = writeRidgeSpectrumFigure(spectrum.frequencyMHz, ...
    {[spectrum.errorPSDdB(:, 1), spectrum.fixedErrorPSDdB(:, 1:3)], ...
    [spectrum.errorPSDdB(:, 2), spectrum.fixedErrorPSDdB(:, 4:6)]}, ...
    false, 'Error PSD relative to target peak (dB)', ridgeErrorFile);

results.selectedPointDirectory = string(directory);
results.outputSpectrumFigure = string(outputFile);
results.errorSpectrumFigure = string(errorFile);
results.ridgeOutputSpectrumFigure = string(ridgeOutputFile);
results.ridgeErrorSpectrumFigure = string(ridgeErrorFile);
results.spectrumConfig = spectrum.config;
results.spectrumConfig.sampleRateSource = string(results.sampleRateSource);
results.spectrumConfig.outputCurveCount = outputCurveCount;
results.spectrumConfig.errorCurveCount = errorCurveCount;
results.spectrumConfig.ridgeOutputPanelCurveCounts = ridgeOutputCounts;
results.spectrumConfig.ridgeErrorPanelCurveCounts = ridgeErrorCounts;
end

function values = fixedPredictionMatrix(predictions)
families = {'complexGMP','pnIQ'};
fields = {'lambda1e3','lambda1e4','lambda1e5'};
values = [];
for familyIndex = 1:numel(families)
    family = families{familyIndex};
    if ~isfield(predictions, family) || ...
            ~all(isfield(predictions.(family), fields))
        error('writeSelectedPointSpectra:MissingFixedPredictions', ...
            'Two fixed-ridge families with three lambdas are required.');
    end
    for fieldIndex = 1:numel(fields)
        values(:, end + 1) = ...
            predictions.(family).(fields{fieldIndex})(:); %#ok<AGROW>
    end
end
end

function curveCount = writeSpectrumFigure( ...
    frequencyMHz, values, labels, yLabel, filename)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 940 600]);
cleanup = onCleanup(@() closeIfValid(figureHandle));
lines = plot(frequencyMHz, values, 'LineWidth', 1.15);
curveCount = numel(lines);
if curveCount ~= numel(labels)
    error('writeSelectedPointSpectra:CurveCountMismatch', ...
        'The plotted curve count does not match the supplied labels.');
end
for index = 1:curveCount
    lines(index).DisplayName = labels(index);
end
grid on;
xlabel('Frequency (MHz)');
ylabel(yLabel);
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', min(curveCount, 4), 'Interpreter', 'none');
exportFigure(figureHandle, filename);
clear cleanup;
end

function curveCounts = writeRidgeSpectrumFigure( ...
    frequencyMHz, panels, includeTarget, yLabel, filename)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1120 780]);
cleanup = onCleanup(@() closeIfValid(figureHandle));
layout = tiledlayout(figureHandle, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
titles = ["Complex GMP-DOMP", "PN-IQ PN-DOMP"];
familyColors = lines(2);
lineStyles = ["-", "--", "-.", ":"];
curveCounts = zeros(1, 2);
for panelIndex = 1:2
    axesHandle = nexttile(layout);
    hold(axesHandle, 'on');
    values = panels{panelIndex};
    if includeTarget
        targetLine = plot(axesHandle, frequencyMHz, values(:, 1), ...
            'Color', [0.15 0.15 0.15], 'LineWidth', 1.2, ...
            'DisplayName', 'Target full signal');
        familyValues = values(:, 2:end);
    else
        targetLine = gobjects(0);
        familyValues = values;
    end
    familyLines = gobjects(1, 4);
    labels = [titles(panelIndex) + " principal", ...
        titles(panelIndex) + ", lambda=1e-3", ...
        titles(panelIndex) + ", lambda=1e-4", ...
        titles(panelIndex) + ", lambda=1e-5"];
    for lineIndex = 1:4
        familyLines(lineIndex) = plot(axesHandle, frequencyMHz, ...
            familyValues(:, lineIndex), ...
            'Color', familyColors(panelIndex, :), ...
            'LineStyle', lineStyles(lineIndex), ...
            'LineWidth', 1 + (lineIndex == 1), ...
            'DisplayName', labels(lineIndex));
    end
    curveCounts(panelIndex) = numel(targetLine) + numel(familyLines);
    grid(axesHandle, 'on');
    xlabel(axesHandle, 'Frequency (MHz)');
    ylabel(axesHandle, yLabel);
    title(axesHandle, titles(panelIndex), 'Interpreter', 'none');
    legend(axesHandle, 'Location', 'eastoutside', 'Interpreter', 'none');
end
exportFigure(figureHandle, filename);
clear cleanup;
end

function exportFigure(figureHandle, filename)
temporary = [tempname(fileparts(filename)) '.png'];
temporaryCleanup = onCleanup(@() deleteIfPresent(temporary));
exportgraphics(figureHandle, temporary, 'Resolution', 160, ...
    'BackgroundColor', 'white');
[moved, message] = movefile(temporary, filename, 'f');
if ~moved
    error('writeSelectedPointSpectra:FigureMoveFailed', ...
        'Could not install %s: %s', filename, message);
end
clear temporaryCleanup;
end

function deleteIfPresent(filename)
if isfile(filename)
    delete(filename);
end
end

function closeIfValid(handle)
if isgraphics(handle)
    close(handle);
end
end
