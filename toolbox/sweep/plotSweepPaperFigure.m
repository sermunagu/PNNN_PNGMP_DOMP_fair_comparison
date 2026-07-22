function files = plotSweepPaperFigure(results, fixed, xVariable, ...
    xLabel, baseFilename, options)
% plotSweepPaperFigure - Render one canonical sweep figure in paper style.

if nargin < 6 || isempty(options)
    options = struct();
end
options = applyDefaults(options);
style = getIEEEPaperStyle();
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 7.0 4.5]);
cleanup = onCleanup(@() close(figureHandle));
axesHandle = axes(figureHandle);
hold(axesHandle, 'on');

models = [options.names.complexGMPDOMP, ...
    options.names.pniqGMP, options.names.pnnn];
labels = models;
mainLines = gobjects(3, 1);
shownY = zeros(0, 1);
for index = 1:3
    rows = string(results.Model) == models(index);
    x = double(results.(xVariable)(rows));
    y = double(results.(options.metricVariable)(rows));
    [x, order] = sort(x);
    y = y(order);
    mainLines(index) = plot(axesHandle, x, y, ...
        'Color', style.mainColors(index, :), ...
        'LineStyle', style.mainLineStyles{index}, ...
        'Marker', style.mainMarkers{index}, ...
        'MarkerIndices', sweepMarkerIndices(numel(x)), ...
        'MarkerSize', style.markerSize, ...
        'LineWidth', style.mainLineWidth, ...
        'DisplayName', labels(index));
    shownY = [shownY; y(:)]; %#ok<AGROW>
end

if options.includeFixed
    fixedModels = [options.names.complexGMPDOMP, options.names.pniqGMP];
    for modelIndex = 1:2
        for lambdaIndex = 1:numel(options.fixedLambdas)
            rows = string(fixed.Model) == fixedModels(modelIndex) & ...
                fixed.FixedLambda == options.fixedLambdas(lambdaIndex);
            x = double(fixed.(xVariable)(rows));
            y = double(fixed.(options.metricVariable)(rows));
            [x, order] = sort(x);
            y = y(order);
            plot(axesHandle, x, y, ...
                'Color', style.mainColors(modelIndex, :), ...
                'LineStyle', style.ridgeLineStyles{lambdaIndex}, ...
                'LineWidth', style.ridgeLineWidth, ...
                'DisplayName', fixedModels(modelIndex) + ...
                    ", \lambda=" + compose('%g', ...
                    options.fixedLambdas(lambdaIndex)));
            shownY = [shownY; y(:)]; %#ok<AGROW>
        end
    end
end

if options.annotateSelected
    selectedParameters = double(options.selection.selectedParameters);
    selectedX = zeros(3, 1);
    selectedY = zeros(3, 1);
    for modelIndex = 1:3
        row = string(results.Model) == models(modelIndex) & ...
            results.ActualRealParameters == selectedParameters;
        if nnz(row) ~= 1
            error('plotSweepPaperFigure:MissingSelectedBudget', ...
                ['The selected common budget must occur once in each ' ...
                'principal family.']);
        end
        selectedX(modelIndex) = double(results.(xVariable)(row));
        selectedY(modelIndex) = ...
            double(results.(options.metricVariable)(row));
    end
    if strcmp(xVariable, 'ActualRealParameters')
        xline(axesHandle, selectedParameters, ':', ...
            'Color', style.selectedRed, 'LineWidth', 1.1, ...
            'HandleVisibility', 'off');
    end
    plot(axesHandle, selectedX, selectedY, 'o', ...
        'LineStyle', 'none', 'Color', style.selectedRed, ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', style.markerSize + 2, 'LineWidth', 0.8, ...
        'HandleVisibility', 'off', ...
        'DisplayName', 'All three models used in selection');
    plot(axesHandle, selectedX(2), selectedY(2), 'o', ...
        'LineStyle', 'none', 'Color', style.selectedRed, ...
        'MarkerFaceColor', style.selectedRed, ...
        'MarkerSize', 3.5, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off', ...
        'DisplayName', options.names.pniqGMP + ...
            " at selected common budget");
end

grid(axesHandle, 'on');
box(axesHandle, 'on');
xlabel(axesHandle, xLabel);
ylabel(axesHandle, options.metricLabel);
set(axesHandle, 'FontName', style.fontName, ...
    'FontSize', style.fontSize, 'LineWidth', style.axisLineWidth);
if options.isNMSE
    finiteY = shownY(isfinite(shownY));
    lowerLimit = floor(min(finiteY)/5)*5;
    if lowerLimit >= -30
        lowerLimit = -35;
    end
    ylim(axesHandle, [lowerLimit -30]);
elseif options.useLogWhenPositive
    finiteY = shownY(isfinite(shownY));
    if all(finiteY > 0)
        set(axesHandle, 'YScale', 'log');
    else
        set(axesHandle, 'YScale', 'linear');
        fprintf(['[Output] MaxAbsRealParameter contains nonpositive values; ' ...
            'using a linear y-axis.\n']);
    end
end
legend(axesHandle, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'NumColumns', 3, ...
    'Interpreter', 'tex', 'FontName', style.fontName, ...
    'FontSize', style.fontSize - 1);
exportOptions = options.exportOptions;
markerStep = max(1, ceil(nnz(string(results.Model) == ...
    options.names.complexGMPDOMP)/17));
exportOptions.tikzExtraAxisOptions = { ...
    'clip=true', ...
    'clip mode=individual', ...
    'legend columns=3', ...
    'legend style={at={(0.5,-0.16)},anchor=north}', ...
    sprintf('every axis plot/.append style={mark repeat=%d}', markerStep)};
files = exportPaperFigure(figureHandle, baseFilename, exportOptions);
clear cleanup;
end

function indices = sweepMarkerIndices(count)
step = max(1, ceil(count/17));
indices = unique([1:step:count, count]);
end

function options = applyDefaults(options)
if ~isfield(options, 'metricVariable')
    options.metricVariable = 'FullSignalNMSEdB';
end
if ~isfield(options, 'metricLabel')
    options.metricLabel = 'Validation NMSE (dB)';
end
if ~isfield(options, 'names')
    error('plotSweepPaperFigure:MissingNames', ...
        'Canonical public names must be supplied in options.names.');
end
if ~isfield(options, 'includeFixed')
    options.includeFixed = false;
end
if ~isfield(options, 'fixedLambdas')
    options.fixedLambdas = [1e-3 1e-4 1e-5];
end
if ~isfield(options, 'isNMSE')
    options.isNMSE = false;
end
if ~isfield(options, 'useLogWhenPositive')
    options.useLogWhenPositive = false;
end
if ~isfield(options, 'annotateSelected')
    options.annotateSelected = false;
end
if ~isfield(options, 'selection')
    options.selection = struct();
end
if ~isfield(options, 'exportOptions')
    options.exportOptions = struct();
end
end
