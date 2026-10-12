function results = run_selected_comparison(selectedParameters, sweep)
% run_selected_comparison - Present one parameter-sweep point.
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
cfg = sweep.config;
addpath(cfg.paper.matlab2tikzSource);
directory = char(sweep.resultDirectory);
fixedLambdas = cfg.fixedRidgeLambdas(:);
automaticSelection = sweep.selection;
modelNames = [cfg.names.complexGMPDOMP; cfg.names.pniqGMP; cfg.names.pnnn];
selectedRows = sweep.results([],:);
for model = modelNames.'
    row = string(sweep.results.Model) == model & ...
        sweep.results.TargetRealParameters == selectedParameters;
    selectedRows = [selectedRows; sweep.results(row,:)]; %#ok<AGROW>
end
stored = load(fullfile(directory, 'linear_sweep.mat'), 'linear');
linear = stored.linear;
stored = load(fullfile(directory, sprintf('pnnn_target_%04d.mat', ...
    selectedParameters)), 'point');
pnnn = stored.point;
complexIndex = find(linear.complexTable.TargetRealParameters == selectedParameters);
pniqIndex = find(linear.pniqTable.TargetRealParameters == selectedParameters);
complexPrediction = linear.predictions.complexFull(:, complexIndex);
pniqPrediction = linear.predictions.pniqFull(:, pniqIndex);
pnnnPrediction = pnnn.fullSignalPrediction(:);

%% Use the preprocessed measurement and alignment saved by this run
stored = load(fullfile(directory, 'experiment.mat'), 'x', 'y', 'split');
x = stored.x;
y = stored.y;
split = stored.split;
fullSignalIndices = double(split.fullSignalIndices(:));
modelInput = x(fullSignalIndices);
targetFullSignal = y(fullSignalIndices);
measurement = load(cfg.measurementFile, 'fs', 'info_signal');
sampleRateHz = resolveMeasurementSampleRate(measurement);

%% Refit fixed Ridge on the saved supports to obtain its predictions
selectedLinear = struct( ...
    'complexTable', linear.complexTable(complexIndex, :), ...
    'pniqTable', linear.pniqTable(pniqIndex, :), ...
    'paths', linear.paths, 'pniqPathMap', linear.pniqPathMap);
selectedConfig = cfg;
selectedConfig.sweep.parameterGrid = selectedParameters;
evaluated = run_fixed_ridge_sweep(x, y, split, selectedConfig, selectedLinear, true);
fixedTable = orderFixedRows(evaluated.table, selectedParameters, fixedLambdas, cfg.names);
fixedPredictionMatrix = [evaluated.predictions.complexFull, evaluated.predictions.pniqFull];

%% Build the public three-row comparison
comparisonTable = selectedRows;
supportCount = selectedParameters/2;
pniqFeatures = linear.paths.pniq(1:supportCount);
completeComparisonTable = buildCompleteComparisonTable( ...
    comparisonTable, fixedTable, cfg.names);
results = struct('selectedParameters', selectedParameters, ...
    'automaticSelectedParameters', automaticSelection.selectedParameters, ...
    'operatingPointSelection', automaticSelection, ...
    'criterionName', automaticSelection.criterionName, ...
    'stabilizationWindowParameters', ...
        automaticSelection.stabilizationWindowParameters, ...
    'stabilizationToleranceDb', ...
        automaticSelection.stabilizationToleranceDb, ...
    'comparisonTable', comparisonTable, ...
    'fixedLambdaComparisonTable', fixedTable, ...
    'completeComparisonTable', completeComparisonTable, ...
    'linearSupports', struct( ...
        'complex', linear.paths.complex(1:supportCount), ...
        'pniqFeatures', pniqFeatures, ...
        'pniqComplex', unique(linear.pniqPathMap.SourceRegressorIndex( ...
            1:supportCount), 'stable')), ...
    'sampleRateHz', sampleRateHz, ...
    'evaluationProtocol', sweep.evaluationProtocol);

%% Compute one shared Welch spectrum and export four paper figures
mainPredictions = [complexPrediction, pniqPrediction, pnnnPrediction];
spectrum = computeSelectedPointSpectra(targetFullSignal, mainPredictions, ...
    sampleRateHz, fixedPredictionMatrix);
selectedDirectory = fullfile(directory, ...
    sprintf('selected_point_%04d', selectedParameters));
if ~isfolder(selectedDirectory)
    mkdir(selectedDirectory);
end
completeTableFile = fullfile(selectedDirectory, ...
    'selected_complete_comparison.csv');
writetable(completeComparisonTable, completeTableFile);
results.completeComparisonTableFile = string(completeTableFile);
exportOptions = struct();
[timeDomainFigureFiles, timeDomainWindow] = ...
    exportSelectedTimeDomainFigures(targetFullSignal, complexPrediction, ...
    pniqPrediction, selectedDirectory, cfg.names, exportOptions);
results.timeDomainFigureFiles = timeDomainFigureFiles;
results.timeDomainRealFigure = timeDomainFigureFiles.real.png;
results.timeDomainImaginaryFigure = timeDomainFigureFiles.imaginary.png;
results.timeDomainWindow = timeDomainWindow;
ampmFigureFiles = exportSelectedAMPMFigure(modelInput, targetFullSignal, ...
    complexPrediction, pniqPrediction, selectedDirectory, cfg.names, ...
    exportOptions);
results.ampmFigureFiles = ampmFigureFiles;
results.ampmFigure = ampmFigureFiles.png;
phaseErrorFigureFiles = exportSelectedPhaseErrorFigure(modelInput, ...
    targetFullSignal, complexPrediction, pniqPrediction, selectedDirectory, ...
    cfg.names, exportOptions);
results.phaseErrorFigureFiles = phaseErrorFigureFiles;
results.phaseErrorFigure = phaseErrorFigureFiles.png;
figureFiles = exportSelectedSpectrumFigures(spectrum, selectedDirectory, ...
    fixedLambdas, cfg.names, exportOptions);
results.spectrumFigureFiles = figureFiles;
results.outputSpectrumFigure = figureFiles.output.png;
results.errorSpectrumFigure = figureFiles.error.png;
results.ridgeOutputSpectrumFigure = figureFiles.ridgeOutput.png;
results.ridgeErrorSpectrumFigure = figureFiles.ridgeError.png;
results.spectrumConfig = spectrum.config;
fprintf('\n=== Complete tutor comparison (principal and fixed Ridge) ===\n');
disp(completeComparisonTable);
fprintf('Complete tutor table: %s\n', completeTableFile);
fprintf('Reused linear and PNNN checkpoints; fixed Ridge refitted on saved supports\n');
fprintf('Selected-point figures: %s\n', selectedDirectory);
end

function files = exportSelectedAMPMFigure(modelInput, target, ...
    complexPrediction, pniqPrediction, directory, names, exportOptions)
% Export the measured and modeled full-signal AM/PM point clouds.

normalizedInputAmplitude = abs(modelInput) / max(abs(modelInput));
targetAMPM = 180/pi * fase_pmpi(angle(target) - angle(modelInput));
complexAMPM = 180/pi * fase_pmpi( ...
    angle(complexPrediction) - angle(modelInput));
pniqAMPM = 180/pi * fase_pmpi(angle(pniqPrediction) - angle(modelInput));
plotIndices = selectDeterministicPlotIndices(numel(modelInput), 5000);

style = getIEEEPaperStyle();
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 7.0 4.2]);
cleanup = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
plot(axesHandle, normalizedInputAmplitude(plotIndices), ...
    targetAMPM(plotIndices), '.', 'Color', 'k', 'MarkerSize', 4);
hold(axesHandle, 'on');
plot(axesHandle, normalizedInputAmplitude(plotIndices), ...
    complexAMPM(plotIndices), '.', 'Color', [0 0.45 0.74], ...
    'MarkerSize', 4);
plot(axesHandle, normalizedInputAmplitude(plotIndices), ...
    pniqAMPM(plotIndices), '.', 'Color', [0.85 0.33 0.10], ...
    'MarkerSize', 4);
ylabel(axesHandle, 'Phase Shift (deg)');
xlabel(axesHandle, 'Normalized Input Amplitude');
title(axesHandle, 'Measured and Modeled AM/PM Characteristics', ...
    'Interpreter', 'none', 'FontWeight', 'normal');
legend(axesHandle, names.measuredOutput, names.complexGMPDOMP, ...
    names.pniqGMP, ...
    'Location', 'northeast', 'FontName', style.fontName, ...
    'FontSize', style.fontSize - 1);
grid(axesHandle, 'on');
box(axesHandle, 'on');
set(axesHandle, 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'LineWidth', style.axisLineWidth);
files = exportPaperFigure(figureHandle, ...
    fullfile(directory, 'selected_ampm'), exportOptions);
clear cleanup;
end

function files = exportSelectedPhaseErrorFigure(modelInput, target, ...
    complexPrediction, pniqPrediction, directory, names, exportOptions)
% Export direct target-minus-prediction phase errors versus input amplitude.

normalizedInputAmplitude = abs(modelInput) / max(abs(modelInput));
complexPhaseError = 180/pi * fase_pmpi( ...
    angle(target) - angle(complexPrediction));
pniqPhaseError = 180/pi * fase_pmpi( ...
    angle(target) - angle(pniqPrediction));
plotIndices = selectDeterministicPlotIndices(numel(modelInput), 5000);

style = getIEEEPaperStyle();
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 7.0 4.2]);
cleanup = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
plot(axesHandle, normalizedInputAmplitude(plotIndices), ...
    complexPhaseError(plotIndices), '.', ...
    'Color', [0 0.45 0.74], 'MarkerSize', 4);
hold(axesHandle, 'on');
plot(axesHandle, normalizedInputAmplitude(plotIndices), ...
    pniqPhaseError(plotIndices), '.', ...
    'Color', [0.85 0.33 0.10], 'MarkerSize', 4);
ylabel(axesHandle, 'Phase error, target - prediction (deg)')
xlabel(axesHandle, 'Normalized Input Amplitude');
title(axesHandle, 'Model Phase Error Characteristics', ...
    'Interpreter', 'none', 'FontWeight', 'normal');
legend(axesHandle, names.complexGMPDOMP, names.pniqGMP, ...
    'Location', 'northeast', 'FontName', style.fontName, ...
    'FontSize', style.fontSize - 1);
grid(axesHandle, 'on');
box(axesHandle, 'on');
set(axesHandle, 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'LineWidth', style.axisLineWidth);
files = exportPaperFigure(figureHandle, ...
    fullfile(directory, 'selected_phase_error'), exportOptions);
clear cleanup;
end

function indices = selectDeterministicPlotIndices(sampleCount, maximumPoints)
% Keep TikZ exports bounded without changing any computed result.

plotCount = min(sampleCount, maximumPoints);
indices = unique(round(linspace(1, sampleCount, plotCount))).';
end

function phase = fase_pmpi(phase)
phase = mod(phase + pi, 2*pi) - pi;
end

function complete = buildCompleteComparisonTable(mainRows, fixedRows, names)
% Combine the three principal rows and six fixed-Ridge references.

familyOrder = [names.complexGMPDOMP, names.pniqGMP, names.pnnn];
complete = table();
for family = familyOrder
    if family == names.complexGMPDOMP
        mainIndex = 1;
    elseif family == names.pniqGMP
        mainIndex = 2;
    else
        mainIndex = 3;
    end
    principal = mainRows(mainIndex, :);
    variant = "Principal / lambda=0";
    if family == names.pnnn
        variant = "Principal";
    end
    principalRow = table(family, variant, false, ...
        double(principal.SelectedLambda), ...
        double(principal.TargetRealParameters), ...
        double(principal.ActualRealParameters), ...
        double(principal.InternalValidationNMSEdB), ...
        double(principal.IdentificationNMSEdB), ...
        double(principal.FullSignalNMSEdB), ...
        double(principal.FLOPsPerSample), ...
        double(principal.ActiveWeights), ...
        double(principal.ActiveBiases), ...
        double(principal.WeightSparsityPercent), ...
        double(principal.FineTuneEpochs), ...
        double(principal.MaxAbsRealParameter), ...
        'VariableNames', {'Family','Variant','UsesFixedRidge','Lambda', ...
        'TargetRealParameters','ActualRealParameters', ...
        'InternalValidationNMSEdB','IdentificationNMSEdB', ...
        'FullSignalNMSEdB','FLOPsPerSample','ActiveWeights', ...
        'ActiveBiases','WeightSparsityPercent','FineTuneEpochs', ...
        'MaxAbsRealParameter'});
    if isempty(complete.Properties.VariableNames)
        complete = principalRow;
    else
        complete = [complete; principalRow]; %#ok<AGROW>
    end
    if family == names.pnnn
        continue;
    end
    familyFixed = fixedRows(string(fixedRows.Model) == family, :);
    for row = 1:height(familyFixed)
        ridgeRow = table(family, "Fixed Ridge", true, ...
            double(familyFixed.FixedLambda(row)), ...
            double(familyFixed.TargetRealParameters(row)), ...
            double(familyFixed.ActualRealParameters(row)), NaN, ...
            double(familyFixed.IdentificationNMSEdB(row)), ...
            double(familyFixed.FullSignalNMSEdB(row)), ...
            double(familyFixed.FLOPsPerSample(row)), ...
            double(familyFixed.ActualRealParameters(row)), 0, NaN, NaN, ...
            double(familyFixed.MaxAbsRealParameter(row)), ...
            'VariableNames', principalRow.Properties.VariableNames);
        complete = [complete; ridgeRow]; %#ok<AGROW>
    end
end
end

function rows = orderFixedRows(value, target, fixedLambdas, names)
models = [names.complexGMPDOMP; names.pniqGMP];
rows = value([],:);
for model = models.'
    for lambda = fixedLambdas.'
        match = string(value.Model) == model & ...
            value.TargetRealParameters == target & value.FixedLambda == lambda;
        rows = [rows; value(match,:)]; %#ok<AGROW>
    end
end
end

function [files, window] = exportSelectedTimeDomainFigures( ...
    target, complexPrediction, pniqPrediction, directory, names, exportOptions)
% Export deterministic central time-domain views for the two linear models.

target = target(:);
complexPrediction = complexPrediction(:);
pniqPrediction = pniqPrediction(:);
sampleCount = numel(target);
windowLength = min(1000, sampleCount);
firstSample = floor((sampleCount - windowLength)/2) + 1;
lastSample = firstSample + windowLength - 1;
indices = (firstSample:lastSample).';
signals = [target(indices), complexPrediction(indices), ...
    pniqPrediction(indices)];

window = struct('firstAlignedSample', firstSample, ...
    'lastAlignedSample', lastSample, 'sampleCount', windowLength);
files = struct();
files.real = exportOneTimeDomainFigure(indices, signals, 'real', ...
    fullfile(directory, 'selected_time_domain_real'), names, exportOptions);
files.imaginary = exportOneTimeDomainFigure(indices, signals, 'imaginary', ...
    fullfile(directory, 'selected_time_domain_imaginary'), names, ...
    exportOptions);
end

function files = exportOneTimeDomainFigure(sampleIndices, signals, component, ...
    baseFilename, names, exportOptions)
style = getIEEEPaperStyle();
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 7.16 4.6]);
cleanup = onCleanup(@() close(figureHandle));
layout = tiledlayout(figureHandle, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
layout.OuterPosition = [0 0.22 1 0.74];

leftAxes = nexttile(layout);
legendLines = plotTimeDomainPanel(leftAxes, sampleIndices, abs(signals), ...
    style);
title(leftAxes, '(a) Magnitude', 'Interpreter', 'none', ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
ylabel(leftAxes, 'Magnitude');

rightAxes = nexttile(layout);
if component == "real"
    componentValues = real(signals);
    componentTitle = '(b) Real part';
    componentLabel = 'Real part';
else
    componentValues = imag(signals);
    componentTitle = '(b) Imaginary part';
    componentLabel = 'Imaginary part';
end
plotTimeDomainPanel(rightAxes, sampleIndices, componentValues, style);
title(rightAxes, componentTitle, 'Interpreter', 'none', ...
    'FontSize', style.fontSize, 'FontWeight', 'normal');
ylabel(rightAxes, componentLabel);

labels = [names.measuredOutput, names.complexGMPDOMP, names.pniqGMP];
legendHandle = legend(leftAxes, legendLines, labels, ...
    'Orientation', 'horizontal', 'NumColumns', 3, ...
    'Interpreter', 'none', 'FontName', style.fontName, ...
    'FontSize', style.fontSize - 1);
legendHandle.Units = 'normalized';
legendHandle.Position = [0.10 0.03 0.80 0.10];
legendHandle.Units = 'pixels';
axesHandles = findobj(figureHandle, 'Type', 'axes');
for axesIndex = 1:numel(axesHandles)
    axesHandles(axesIndex).Units = 'pixels';
end

markerStep = max(1, floor(numel(sampleIndices)/14));
exportOptions.tikzExtraAxisOptions = { ...
    sprintf('every axis plot/.append style={mark repeat=%d}', markerStep)};
files = exportPaperFigure(figureHandle, baseFilename, exportOptions);
clear cleanup;
end

function lines = plotTimeDomainPanel(axesHandle, sampleIndices, values, style)
hold(axesHandle, 'on');
markerIndices = unique([1:max(1, floor(numel(sampleIndices)/14)): ...
    numel(sampleIndices), numel(sampleIndices)]);
lines = gobjects(3, 1);
lines(1) = plot(axesHandle, sampleIndices, values(:, 1), ...
    'Color', style.targetGray, 'LineStyle', '-', ...
    'LineWidth', style.targetLineWidth);
for familyIndex = 1:2
    lines(familyIndex + 1) = plot(axesHandle, sampleIndices, ...
        values(:, familyIndex + 1), ...
        'Color', style.mainColors(familyIndex, :), ...
        'LineStyle', style.mainLineStyles{familyIndex}, ...
        'Marker', style.mainMarkers{familyIndex}, ...
        'MarkerIndices', markerIndices, ...
        'MarkerSize', style.markerSize - 1, ...
        'LineWidth', style.mainLineWidth);
end
grid(axesHandle, 'on');
box(axesHandle, 'on');
xlabel(axesHandle, 'Aligned full-signal sample');
set(axesHandle, 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'LineWidth', style.axisLineWidth);
end
