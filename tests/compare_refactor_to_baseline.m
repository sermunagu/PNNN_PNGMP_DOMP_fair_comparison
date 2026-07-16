% Script: compare_refactor_to_baseline
% Compare the newest complete run with the frozen 20260715_142905 baseline.
% Wall-clock times and corrected descriptive prose are intentionally excluded.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
resultsRoot = fullfile(projectRoot, 'results', ...
    'full_signal_domp_comparison');
baselineDirectory = fullfile(resultsRoot, '20260715_142905');
entries = dir(resultsRoot);
entries = entries([entries.isdir] & ...
    ~ismember({entries.name}, {'.','..','20260715_142905'}));
entries = entries(arrayfun(@(entry) isfile(fullfile(entry.folder, ...
    entry.name, 'full_signal_predictions.mat')), entries));
[~, order] = sort([entries.datenum], 'descend');
assert(~isempty(order), 'No complete candidate run was found.');
candidateDirectory = fullfile(entries(order(1)).folder, ...
    entries(order(1)).name);

predictionRelativeTolerance = 1e-12;
numericTolerance = 1e-11;

baselineSplit = load(fullfile(baselineDirectory, 'split_indices.mat'));
candidateSplit = load(fullfile(candidateDirectory, 'split_indices.mat'));
splitFields = ["identification_indices", "full_signal_indices", ...
    "internal_train_indices", "internal_validation_indices"];
for fieldName = splitFields
    assert(isequal(baselineSplit.(fieldName), candidateSplit.(fieldName)), ...
        'Split mismatch in %s.', fieldName);
end
assert(isequaln(baselineSplit.split, candidateSplit.split), ...
    'The stored split protocol differs from the baseline.');

baselineDOMP = load(fullfile(baselineDirectory, 'domp_supports.mat'));
candidateDOMP = load(fullfile(candidateDirectory, 'domp_supports.mat'));
supportFields = string(fieldnames(baselineDOMP.supports));
assert(isequal(supportFields, string(fieldnames(candidateDOMP.supports))), ...
    'The set of stored DOMP supports differs from the baseline.');
for fieldName = supportFields.'
    assert(isequal(baselineDOMP.supports.(fieldName), ...
        candidateDOMP.supports.(fieldName)), ...
        'Ordered DOMP support mismatch in %s.', fieldName);
end

baselineResults = readtable(fullfile(baselineDirectory, ...
    'comparison_results.csv'), 'TextType', 'string');
candidateResults = readtable(fullfile(candidateDirectory, ...
    'comparison_results.csv'), 'TextType', 'string');
candidateHistoricalResults = historicalRows(candidateResults);
compareTables(baselineResults, candidateHistoricalResults, ...
    ["TrainingTimeSeconds", "AdditionalOperationsPerSample"], ...
    numericTolerance);

exactTables = ["selected_hyperparameters.csv", "parameter_summary.csv", ...
    "complexity_flops.csv", "flop_convention.csv"];
for filename = exactTables
    baselineTable = readtable(fullfile(baselineDirectory, filename), ...
        'TextType', 'string');
    candidateTable = readtable(fullfile(candidateDirectory, filename), ...
        'TextType', 'string');
    baselineTable = normalizeHistoricalActivationName(baselineTable);
    candidateTable = normalizeHistoricalActivationName(candidateTable);
    if ismember('RegularizationMode', candidateTable.Properties.VariableNames)
        candidateTable = historicalRows(candidateTable);
    end
    compareTables(baselineTable, candidateTable, ...
        strings(0, 1), numericTolerance);
end

baselineMetadata = load(fullfile(baselineDirectory, ...
    'comparison_results.mat'), 'comparison_metadata');
candidateMetadata = load(fullfile(candidateDirectory, ...
    'comparison_results.mat'), 'comparison_metadata');
baselineSparse = baselineMetadata.comparison_metadata.n12Sparse;
candidateSparse = candidateMetadata.comparison_metadata.n12Sparse;
assert(isequaln(baselineSparse.pruningStats, candidateSparse.pruningStats), ...
    'Sparse pruning and mask-integrity summaries differ.');
assert(isequaln(baselineSparse.maskIntegrityAfterPruning, ...
    candidateSparse.maskIntegrityAfterPruning) && ...
    isequaln(baselineSparse.maskIntegrityAfterFineTune, ...
    candidateSparse.maskIntegrityAfterFineTune), ...
    'Sparse mask-integrity checks differ.');

baselinePredictions = load(fullfile(baselineDirectory, ...
    'full_signal_predictions.mat'));
candidatePredictions = load(fullfile(candidateDirectory, ...
    'full_signal_predictions.mat'));
assert(isequal(baselinePredictions.target_identification, ...
    candidatePredictions.target_identification), ...
    'Identification targets differ from the baseline.');
assert(isequal(baselinePredictions.target_full_signal, ...
    candidatePredictions.target_full_signal), ...
    'Full-signal targets differ from the baseline.');

modelFields = string(fieldnames(baselinePredictions.full_signal_predictions));
assert(isequal(modelFields, ...
    string(fieldnames(candidatePredictions.full_signal_predictions))), ...
    'The set of stored prediction models differs from the baseline.');

summary = table('Size', [numel(modelFields), 5], ...
    'VariableTypes', ["string", repmat("double", 1, 4)], ...
    'VariableNames', ["ModelField", "IdentificationRelativeError", ...
    "FullSignalRelativeError", "IdentificationMaximumAbsoluteError", ...
    "FullSignalMaximumAbsoluteError"]);

for modelIndex = 1:numel(modelFields)
    fieldName = modelFields(modelIndex);
    candidateIdentification = ...
        candidatePredictions.identification_predictions.(fieldName);
    candidateFull = candidatePredictions.full_signal_predictions.(fieldName);
    if isstruct(candidateIdentification)
        candidateIdentification = candidateIdentification.selectedRidge;
        candidateFull = candidateFull.selectedRidge;
    end
    [identificationRelative, identificationMaximum] = predictionError( ...
        baselinePredictions.identification_predictions.(fieldName), ...
        candidateIdentification);
    [fullRelative, fullMaximum] = predictionError( ...
        baselinePredictions.full_signal_predictions.(fieldName), ...
        candidateFull);

    summary(modelIndex, :) = {fieldName, identificationRelative, fullRelative, ...
        identificationMaximum, fullMaximum};
end

assert(all(summary.IdentificationRelativeError <= predictionRelativeTolerance), ...
    'At least one identification prediction differs from the baseline.');
assert(all(summary.FullSignalRelativeError <= predictionRelativeTolerance), ...
    'At least one full-signal prediction differs from the baseline.');

fprintf('\nRefactor-to-baseline comparison: PASS\n');
fprintf('Baseline:  %s\n', baselineDirectory);
fprintf('Candidate: %s\n\n', candidateDirectory);
disp(summary);

function historical = historicalRows(results)
linearRow = results.RegularizationMode == "Validation-selected Ridge";
pnnnRow = results.RegularizationMode == "Not applicable";
historical = removevars(results(linearRow | pnnnRow, :), ...
    'RegularizationMode');
end

function compareTables(baseline, candidate, ignoredVariables, tolerance)
assert(isequal(baseline.Properties.VariableNames, ...
    candidate.Properties.VariableNames), 'Table schemas differ.');
assert(height(baseline) == height(candidate), 'Table heights differ.');

variables = setdiff(string(baseline.Properties.VariableNames), ...
    ignoredVariables, 'stable');
for variableName = variables
    baselineValues = baseline.(variableName);
    candidateValues = candidate.(variableName);
    if isnumeric(baselineValues)
        scale = max(1, max(abs(baselineValues), [], 'all', 'omitnan'));
        difference = abs(baselineValues - candidateValues);
        difference(isnan(baselineValues) & isnan(candidateValues)) = 0;
        assert(all(difference <= tolerance * scale, 'all'), ...
            'Numeric table mismatch in %s.', variableName);
    else
        assert(isequaln(baselineValues, candidateValues), ...
            'Table mismatch in %s.', variableName);
    end
end
end

function value = normalizeHistoricalActivationName(value)
oldName = 'NumELUPerSample';
newName = 'NumActivationEvaluationsPerSample';
if ismember(oldName, value.Properties.VariableNames) && ...
        ~ismember(newName, value.Properties.VariableNames)
    value = renamevars(value, oldName, newName);
end
end

function [relativeError, maximumAbsoluteError] = predictionError(reference, value)
difference = value - reference;
relativeError = norm(difference) / max(norm(reference), eps);
maximumAbsoluteError = max(abs(difference), [], 'all');
end
