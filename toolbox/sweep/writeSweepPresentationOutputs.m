function results = writeSweepPresentationOutputs( ...
    linear, pnnnRows, parameterGrid, directory)
% writeSweepPresentationOutputs - Build the public three-family sweep output.
% This function validates and writes the canonical CSV and comparison figures.
% All targets are handled uniformly under one shared presentation contract.

requiredLinear = {'complexTable','pnTable'};
if ~isstruct(linear) || ~all(isfield(linear, requiredLinear)) || ...
        ~istable(pnnnRows)
    error('writeSweepPresentationOutputs:InvalidInput', ...
        'Linear sweep tables and the sparse PNNN table are required.');
end
targets = double(parameterGrid(:));
if isempty(targets) || any(~isfinite(targets)) || ...
        numel(unique(targets)) ~= numel(targets)
    error('writeSweepPresentationOutputs:InvalidTargets', ...
        'The parameter grid must contain unique finite targets.');
end
results = [linear.complexTable; linear.pnTable; pnnnRows];
validatePresentedResults(results, targets);

directory = string(directory);
if ~isscalar(directory) || ismissing(directory) || strlength(directory) == 0
    error('writeSweepPresentationOutputs:InvalidDirectory', ...
        'The result directory must be a nonempty text scalar.');
end
directoryPath = char(directory);
if ~isfolder(directoryPath)
    mkdir(directoryPath);
end
csvFile = fullfile(directoryPath, 'complexity_sweep.csv');
temporaryCSV = [tempname(directoryPath) '.csv'];
csvCleanup = onCleanup(@() deleteIfPresent(temporaryCSV));
fprintf('[Output] Writing complexity_sweep.csv...\n');
writetable(results, temporaryCSV);
[moved, message] = movefile(temporaryCSV, csvFile, 'f');
if ~moved
    error('writeSweepPresentationOutputs:CSVMoveFailed', ...
        'Could not install the sweep CSV: %s', message);
end
clear csvCleanup;

fprintf('[Output] Writing NMSE-vs-parameters figure...\n');
writeComparisonFigure(results, directoryPath, 'ActualRealParameters', ...
    'Active real parameters', 'comparison_nmse_parameters_sweep.png');
fprintf('[Output] Writing NMSE-vs-FLOPs figure...\n');
writeComparisonFigure(results, directoryPath, 'FLOPsPerSample', ...
    'FLOPs per sample', 'comparison_nmse_flops_sweep.png');
end

function validatePresentedResults(results, targets)
required = {'Model','SweepRole','TargetRealParameters', ...
    'ActualRealParameters','FullSignalNMSEdB','FLOPsPerSample', ...
    'ActiveWeights','ActiveBiases'};
if ~istable(results) || ~all(ismember(required, ...
        results.Properties.VariableNames))
    error('writeSweepPresentationOutputs:InvalidResults', ...
        'The canonical sweep columns are required.');
end
models = ["Complex GMP DOMP sweep", ...
    "Independent PN-IQ PN-DOMP sweep", "Sparse PNNN N12"];
if height(results) ~= 3*numel(targets) || ...
        any(string(results.SweepRole) == "Historical reference") || ...
        any(contains(string(results.Model), "Historical", ...
        'IgnoreCase', true)) || any(~ismember(string(results.Model), models))
    error('writeSweepPresentationOutputs:UnexpectedRows', ...
        'Only the three comparable sweep families may be presented.');
end
expectedModels = sort(models(:));
for index = 1:numel(targets)
    rows = results.TargetRealParameters == targets(index);
    actualModels = sort(string(results.Model(rows)));
    if nnz(rows) ~= 3 || ~isequal(actualModels(:), expectedModels) || ...
            any(results.ActualRealParameters(rows) ~= targets(index))
        error('writeSweepPresentationOutputs:TargetMismatch', ...
            'Each target must contain all three exact-parameter families.');
    end
end
pnnnRows = string(results.Model) == "Sparse PNNN N12";
if any(results.ActiveWeights(pnnnRows) + ...
        results.ActiveBiases(pnnnRows) ~= ...
        results.ActualRealParameters(pnnnRows))
    error('writeSweepPresentationOutputs:PNNNParameterMismatch', ...
        'Sparse PNNN active weights and biases must equal its parameter count.');
end
end

function writeComparisonFigure(results, directory, xVariable, xLabel, filename)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 900 600]);
figureCleanup = onCleanup(@() closeIfValid(figureHandle));
hold on;
models = ["Complex GMP DOMP sweep", ...
    "Independent PN-IQ PN-DOMP sweep", "Sparse PNNN N12"];
labels = ["Complex GMP-DOMP", "PN-IQ PN-DOMP", "Sparse PNNN N12"];
for index = 1:numel(models)
    rows = string(results.Model) == models(index);
    x = results.(xVariable)(rows);
    y = results.FullSignalNMSEdB(rows);
    [x, order] = sort(x);
    plot(x, y(order), '-o', 'DisplayName', labels(index));
end
grid on;
xlabel(xLabel);
ylabel('Full-signal NMSE (dB)');
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Interpreter', 'none');
exportgraphics(figureHandle, fullfile(directory, filename), ...
    'Resolution', 160);
clear figureCleanup;
end

function deleteIfPresent(filename)
if isfile(filename)
    delete(filename);
end
end

function closeIfValid(figureHandle)
if isgraphics(figureHandle)
    close(figureHandle);
end
end
