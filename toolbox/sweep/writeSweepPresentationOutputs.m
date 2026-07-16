function [results, details] = writeSweepPresentationOutputs( ...
    linear, fixedLinear, pnnnRows, parameterGrid, directory)
% writeSweepPresentationOutputs - Build the public three-family sweep output.
% Fixed-ridge linear references remain supplementary to the canonical table.
% All targets are validated before the CSV files and figures are replaced.

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
if ~isstruct(fixedLinear) || ~isfield(fixedLinear, 'table')
    error('writeSweepPresentationOutputs:InvalidFixedLinear', ...
        'The supplementary fixed-ridge table is required.');
end
fixedResults = fixedLinear.table;
validateFixedResults(fixedResults, results, targets);

directory = string(directory);
if ~isscalar(directory) || ismissing(directory) || strlength(directory) == 0
    error('writeSweepPresentationOutputs:InvalidDirectory', ...
        'The result directory must be a nonempty text scalar.');
end
directoryPath = char(directory);
if ~isfolder(directoryPath)
    mkdir(directoryPath);
end
fprintf('[Output] Writing complexity_sweep.csv...\n');
writeTableAtomically(results, directoryPath, 'complexity_sweep.csv');
fprintf('[Output] Writing fixed_lambda_linear_sweep.csv...\n');
writeTableAtomically(fixedResults, directoryPath, ...
    'fixed_lambda_linear_sweep.csv');

fprintf('[Output] Writing NMSE-vs-parameters figure...\n');
details.parameterCurveCount = writeComparisonFigure( ...
    results, fixedResults, directoryPath, 'ActualRealParameters', ...
    'Active real parameters', 'comparison_nmse_parameters_sweep.png', true);
fprintf('[Output] Writing NMSE-vs-FLOPs figure...\n');
details.flopsCurveCount = writeComparisonFigure( ...
    results, table(), directoryPath, 'FLOPsPerSample', ...
    'FLOPs per sample', 'comparison_nmse_flops_sweep.png', false);
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

function validateFixedResults(fixed, principal, targets)
required = {'Model','TargetRealParameters','ActualRealParameters', ...
    'FixedLambda','IdentificationNMSEdB','FullSignalNMSEdB', ...
    'FLOPsPerSample'};
models = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"];
lambdas = [1e-3; 1e-4; 1e-5];
if ~istable(fixed) || ~all(ismember(required, ...
        fixed.Properties.VariableNames)) || ...
        height(fixed) ~= 6*numel(targets) || ...
        ~isequal(sort(unique(string(fixed.Model))), sort(models))
    error('writeSweepPresentationOutputs:InvalidFixedRows', ...
        'The fixed-ridge table must contain two families and three lambdas.');
end
principalModels = ["Complex GMP DOMP sweep"; ...
    "Independent PN-IQ PN-DOMP sweep"];
for modelIndex = 1:numel(models)
    for lambdaIndex = 1:numel(lambdas)
        rows = string(fixed.Model) == models(modelIndex) & ...
            fixed.FixedLambda == lambdas(lambdaIndex);
        if nnz(rows) ~= numel(targets) || ...
                ~isequal(sort(fixed.TargetRealParameters(rows)), ...
                sort(targets))
            error('writeSweepPresentationOutputs:MissingFixedCurve', ...
                'Every fixed family/lambda combination needs all targets.');
        end
    end
    for targetIndex = 1:numel(targets)
        auxiliaryRows = string(fixed.Model) == models(modelIndex) & ...
            fixed.TargetRealParameters == targets(targetIndex);
        principalRow = string(principal.Model) == ...
            principalModels(modelIndex) & ...
            principal.TargetRealParameters == targets(targetIndex);
        if nnz(auxiliaryRows) ~= numel(lambdas) || nnz(principalRow) ~= 1 || ...
                any(fixed.ActualRealParameters(auxiliaryRows) ~= ...
                principal.ActualRealParameters(principalRow)) || ...
                any(fixed.FLOPsPerSample(auxiliaryRows) ~= ...
                principal.FLOPsPerSample(principalRow))
            error('writeSweepPresentationOutputs:FixedComplexityMismatch', ...
                'Fixed ridge must preserve principal support complexity.');
        end
    end
end
end

function writeTableAtomically(value, directory, filename)
finalFile = fullfile(directory, filename);
temporaryFile = [tempname(directory) '.csv'];
cleanup = onCleanup(@() deleteIfPresent(temporaryFile));
writetable(value, temporaryFile);
[moved, message] = movefile(temporaryFile, finalFile, 'f');
if ~moved
    error('writeSweepPresentationOutputs:CSVMoveFailed', ...
        'Could not install %s: %s', filename, message);
end
clear cleanup;
end

function curveCount = writeComparisonFigure( ...
    results, fixed, directory, xVariable, xLabel, filename, includeFixed)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 900 600]);
figureCleanup = onCleanup(@() closeIfValid(figureHandle));
hold on;
models = ["Complex GMP DOMP sweep", ...
    "Independent PN-IQ PN-DOMP sweep", "Sparse PNNN N12"];
labels = ["Complex GMP-DOMP", "PN-IQ PN-DOMP", "Sparse PNNN N12"];
mainHandles = gobjects(numel(models), 1);
for index = 1:numel(models)
    rows = string(results.Model) == models(index);
    x = results.(xVariable)(rows);
    y = results.FullSignalNMSEdB(rows);
    [x, order] = sort(x);
    if includeFixed
        mainHandles(index) = plot(x, y(order), '-o', 'LineWidth', 2, ...
            'DisplayName', labels(index));
    else
        mainHandles(index) = plot(x, y(order), '-o', ...
            'DisplayName', labels(index));
    end
end
curveCount = numel(models);
if includeFixed
    auxiliaryModels = ["Complex GMP-DOMP", "PN-IQ PN-DOMP"];
    lambdas = [1e-3, 1e-4, 1e-5];
    lineStyles = ["--", ":", "-."];
    for modelIndex = 1:numel(auxiliaryModels)
        for lambdaIndex = 1:numel(lambdas)
            rows = string(fixed.Model) == auxiliaryModels(modelIndex) & ...
                fixed.FixedLambda == lambdas(lambdaIndex);
            x = fixed.(xVariable)(rows);
            y = fixed.FullSignalNMSEdB(rows);
            [x, order] = sort(x);
            label = auxiliaryModels(modelIndex) + ", lambda=" + ...
                compose('%g', lambdas(lambdaIndex));
            plot(x, y(order), 'LineStyle', lineStyles(lambdaIndex), ...
                'Color', mainHandles(modelIndex).Color, 'LineWidth', 1, ...
                'DisplayName', label);
            curveCount = curveCount + 1;
        end
    end
end
grid on;
xlabel(xLabel);
ylabel('Full-signal NMSE (dB)');
if includeFixed
    legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
        'NumColumns', 3, 'Interpreter', 'none');
else
    legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Interpreter', 'none');
end
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
