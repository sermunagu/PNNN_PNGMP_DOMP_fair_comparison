function selection = selectOperatingPoint(results, selectionConfig)
% selectOperatingPoint - Select the first jointly stabilized common budget.
% A budget is jointly stabilized when none of the three principal families
% can improve by more than the configured NMSE tolerance within the next
% configured number of real parameters. The lowest such budget is selected.

if nargin < 2 || isempty(selectionConfig)
    selectionConfig = struct();
end
if ~isfield(selectionConfig, 'stabilizationWindowParameters')
    selectionConfig.stabilizationWindowParameters = 100;
end
if ~isfield(selectionConfig, 'stabilizationToleranceDb')
    selectionConfig.stabilizationToleranceDb = 0.20;
end
if ~isfield(selectionConfig, 'sensitivityWindowsParameters')
    selectionConfig.sensitivityWindowsParameters = [80 100 120];
end
if ~isfield(selectionConfig, 'sensitivityTolerancesDb')
    selectionConfig.sensitivityTolerancesDb = [0.15 0.20 0.25];
end

windowParameters = double(selectionConfig.stabilizationWindowParameters);
toleranceDb = double(selectionConfig.stabilizationToleranceDb);
if ~(isscalar(windowParameters) && isfinite(windowParameters) && ...
        windowParameters > 0)
    error('selectOperatingPoint:InvalidWindow', ...
        'The stabilization window must be one positive finite scalar.');
end
if ~(isscalar(toleranceDb) && isfinite(toleranceDb) && toleranceDb >= 0)
    error('selectOperatingPoint:InvalidTolerance', ...
        'The stabilization tolerance must be one nonnegative finite scalar.');
end

criterionName = "joint stabilization minimum-complexity criterion";
modelNames = ["Complex GMP DOMP sweep", ...
    "Independent PN-IQ PN-DOMP sweep", "Sparse PNNN N12"];
required = ["Model", "ActualRealParameters", "FullSignalNMSEdB", ...
    "FLOPsPerSample"];
if ~istable(results) || ~all(ismember(required, ...
        string(results.Properties.VariableNames)))
    error('selectOperatingPoint:InvalidResults', ...
        'The canonical results table is missing required variables.');
end

[budgets, nmse, flops] = collectPrincipalFamilies(results, modelNames);
if any(~isfinite([budgets; nmse(:); flops(:)]))
    error('selectOperatingPoint:NonfiniteData', ...
        'Selection metrics must be finite.');
end
if any(diff(flops, 1, 1) < 0, 'all')
    error('selectOperatingPoint:NonmonotonicFLOPs', ...
        ['FLOPs must be nondecreasing with parameter budget in each ' ...
        'principal family.']);
end

[futureGainDb, hasFullWindow] = calculateFutureGain( ...
    budgets, nmse, windowParameters);
jointlyStabilized = hasFullWindow & ...
    all(futureGainDb <= toleranceDb + 1e-12, 2);
selectedRow = find(jointlyStabilized, 1, 'first');
if isempty(selectedRow)
    error('selectOperatingPoint:NoJointlyStabilizedPoint', ...
        ['No common budget jointly stabilizes all three principal ' ...
        'families for the configured window and tolerance.']);
end
selected = false(numel(budgets), 1);
selected(selectedRow) = true;

WindowUpperParameters = budgets + windowParameters;
WindowUpperParameters(~hasFullWindow) = NaN;
ComplexGMPFutureGainDb = futureGainDb(:, 1);
PNIQFutureGainDb = futureGainDb(:, 2);
SparsePNNNFutureGainDb = futureGainDb(:, 3);
WorstFutureGainDb = max(futureGainDb, [], 2);
WorstFutureGainDb(~hasFullWindow) = NaN;
ComplexGMPFLOPs = flops(:, 1);
PNIQFLOPs = flops(:, 2);
SparsePNNNFLOPs = flops(:, 3);
diagnostics = table(budgets, WindowUpperParameters, hasFullWindow, ...
    ComplexGMPFutureGainDb, PNIQFutureGainDb, SparsePNNNFutureGainDb, ...
    WorstFutureGainDb, jointlyStabilized, selected, ComplexGMPFLOPs, ...
    PNIQFLOPs, SparsePNNNFLOPs, 'VariableNames', { ...
    'ActualRealParameters', 'WindowUpperParameters', 'HasFullWindow', ...
    'ComplexGMPFutureGainDb', 'PNIQFutureGainDb', ...
    'SparsePNNNFutureGainDb', 'WorstFutureGainDb', ...
    'JointlyStabilized', 'Selected', 'ComplexGMPFLOPs', ...
    'PNIQFLOPs', 'SparsePNNNFLOPs'});

sensitivity = buildSensitivityTable(selectionConfig, budgets, nmse);
summarySentence = compose([ ...
    'Using a %d-parameter forward window and a %.2f dB stabilization ' ...
    'tolerance, %d active real parameters is the lowest common budget ' ...
    'for which none of the three principal techniques can improve by ' ...
    'more than %.2f dB within the next %d parameters. The remaining ' ...
    'improvements are %.4f dB for Complex GMP-DOMP, %.4f dB for PN-IQ ' ...
    'PN-DOMP, and %.4f dB for Sparse PNNN N12. Because FLOPs increase ' ...
    'monotonically with parameter budget in all three sweeps, this is ' ...
    'also the minimum-FLOP jointly stabilized common configuration.'], ...
    windowParameters, toleranceDb, budgets(selectedRow), toleranceDb, ...
    windowParameters, futureGainDb(selectedRow, 1), ...
    futureGainDb(selectedRow, 2), futureGainDb(selectedRow, 3));

selection = struct( ...
    'selectedParameters', budgets(selectedRow), ...
    'criterionName', criterionName, ...
    'stabilizationWindowParameters', windowParameters, ...
    'stabilizationToleranceDb', toleranceDb, ...
    'selectedWindowUpperParameters', budgets(selectedRow) + ...
        windowParameters, ...
    'selectedComplexGMPFutureGainDb', futureGainDb(selectedRow, 1), ...
    'selectedPNIQFutureGainDb', futureGainDb(selectedRow, 2), ...
    'selectedSparsePNNNFutureGainDb', futureGainDb(selectedRow, 3), ...
    'selectedWorstFutureGainDb', max(futureGainDb(selectedRow, :)), ...
    'selectedComplexGMPFLOPs', flops(selectedRow, 1), ...
    'selectedPNIQFLOPs', flops(selectedRow, 2), ...
    'selectedSparsePNNNFLOPs', flops(selectedRow, 3), ...
    'diagnosticsTable', diagnostics, ...
    'sensitivityTable', sensitivity, ...
    'summarySentence', string(summarySentence));
end

function [budgets, nmse, flops] = collectPrincipalFamilies(results, modelNames)
budgets = [];
nmse = [];
flops = [];
for modelIndex = 1:numel(modelNames)
    rows = string(results.Model) == modelNames(modelIndex);
    family = results(rows, :);
    if isempty(family)
        error('selectOperatingPoint:MissingFamily', ...
            'Missing principal family: %s.', modelNames(modelIndex));
    end
    [familyBudgets, order] = sort(double(family.ActualRealParameters));
    if numel(unique(familyBudgets)) ~= numel(familyBudgets)
        error('selectOperatingPoint:DuplicateBudget', ...
            'Each family must contain one row per real-parameter budget.');
    end
    if modelIndex == 1
        budgets = familyBudgets;
        nmse = zeros(numel(budgets), numel(modelNames));
        flops = zeros(numel(budgets), numel(modelNames));
    elseif ~isequal(familyBudgets, budgets)
        error('selectOperatingPoint:BudgetMismatch', ...
            'The three principal families must share the same budgets.');
    end
    nmse(:, modelIndex) = double(family.FullSignalNMSEdB(order));
    flops(:, modelIndex) = double(family.FLOPsPerSample(order));
end
end

function [futureGainDb, hasFullWindow] = calculateFutureGain( ...
    budgets, nmse, windowParameters)
rowCount = numel(budgets);
futureGainDb = nan(rowCount, size(nmse, 2));
hasFullWindow = budgets + windowParameters <= budgets(end) + 1e-12;
for row = find(hasFullWindow).'
    futureRows = budgets >= budgets(row) & ...
        budgets <= budgets(row) + windowParameters + 1e-12;
    futureGainDb(row, :) = nmse(row, :) - min(nmse(futureRows, :), [], 1);
end
end

function sensitivity = buildSensitivityTable(config, budgets, nmse)
windows = double(config.sensitivityWindowsParameters(:));
tolerances = double(config.sensitivityTolerancesDb(:));
rowCount = numel(windows)*numel(tolerances);
StabilizationWindowParameters = zeros(rowCount, 1);
StabilizationToleranceDb = zeros(rowCount, 1);
HasJointlyStabilizedPoint = false(rowCount, 1);
SelectedParameters = nan(rowCount, 1);
SelectedComplexGMPFutureGainDb = nan(rowCount, 1);
SelectedPNIQFutureGainDb = nan(rowCount, 1);
SelectedSparsePNNNFutureGainDb = nan(rowCount, 1);
SelectedWorstFutureGainDb = nan(rowCount, 1);
outputRow = 0;
for windowIndex = 1:numel(windows)
    if ~(isfinite(windows(windowIndex)) && windows(windowIndex) > 0)
        error('selectOperatingPoint:InvalidSensitivityWindow', ...
            'Sensitivity windows must be positive and finite.');
    end
    [gainDb, hasFullWindow] = calculateFutureGain( ...
        budgets, nmse, windows(windowIndex));
    for toleranceIndex = 1:numel(tolerances)
        if ~(isfinite(tolerances(toleranceIndex)) && ...
                tolerances(toleranceIndex) >= 0)
            error('selectOperatingPoint:InvalidSensitivityTolerance', ...
                'Sensitivity tolerances must be nonnegative and finite.');
        end
        outputRow = outputRow + 1;
        StabilizationWindowParameters(outputRow) = windows(windowIndex);
        StabilizationToleranceDb(outputRow) = tolerances(toleranceIndex);
        admissible = hasFullWindow & ...
            all(gainDb <= tolerances(toleranceIndex) + 1e-12, 2);
        selectedRow = find(admissible, 1, 'first');
        if isempty(selectedRow)
            continue;
        end
        HasJointlyStabilizedPoint(outputRow) = true;
        SelectedParameters(outputRow) = budgets(selectedRow);
        SelectedComplexGMPFutureGainDb(outputRow) = gainDb(selectedRow, 1);
        SelectedPNIQFutureGainDb(outputRow) = gainDb(selectedRow, 2);
        SelectedSparsePNNNFutureGainDb(outputRow) = gainDb(selectedRow, 3);
        SelectedWorstFutureGainDb(outputRow) = max(gainDb(selectedRow, :));
    end
end
sensitivity = table(StabilizationWindowParameters, ...
    StabilizationToleranceDb, HasJointlyStabilizedPoint, ...
    SelectedParameters, SelectedComplexGMPFutureGainDb, ...
    SelectedPNIQFutureGainDb, SelectedSparsePNNNFutureGainDb, ...
    SelectedWorstFutureGainDb);
end
