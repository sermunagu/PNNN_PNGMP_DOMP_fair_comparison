% Script: run_fair_PNNN_vs_PNGMP_DOMP
% Run the shared 10% identification/full-signal comparison from end to end.
% Linear models run first; PNNN reuses their active-parameter target.

clearvars;
clc;

project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'toolbox', 'comparison'));
addpath(fullfile(project_root, 'toolbox', 'metrics'));
addpath(fullfile(project_root, 'toolbox', 'complexity'));
addpath(fullfile(project_root, 'toolbox', 'domp'));
addpath(fullfile(project_root, 'toolbox', 'splits'));
addpath(fullfile(project_root, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(project_root, 'toolbox', 'pnnn'));
addpath(fullfile(project_root, 'toolbox', 'pnnn', 'pruning'));
addpath(fullfile(project_root, 'toolbox', 'reporting'));

cfg = getFairDOMPComparisonConfig(project_root);
if cfg.warmStart.enabled || cfg.warmStart.useLatestDeploy
    error('run_fair_PNNN_vs_PNGMP_DOMP:WarmStartForbidden', ...
        'The fair PNNN comparison must initialize final networks from scratch.');
end

[data, split, result_directory] = prepareComparison(cfg);
diary(fullfile(result_directory, 'run_log.txt'));
diary_cleanup = onCleanup(@() diary('off'));

fprintf('\n=== Full-signal PNNN versus PN-GMP DOMP comparison ===\n');
fprintf('Result directory: %s\n', result_directory);
fprintf('Mapping: %s (local modeled-block X/Y convention)\n', ...
    cfg.mappingMode);
fprintf(['Internal train=%d | internal validation=%d | ' ...
    'identification=%d | full signal=%d\n'], ...
    numel(split.internalTrainIndices), ...
    numel(split.internalValidationIndices), ...
    numel(split.identificationIndices), ...
    numel(split.fullSignalIndices));
fprintf('Identification is contained in full signal: YES\n');

fprintf('\nFitting six linear DOMP models under the corrected protocol...\n');
linear_study = runPNGMPDOMPStudy(data.x, data.y, split, cfg);
independent_row = ...
    linear_study.comparisonResults.Model == "Independent PN-IQ full";
matched_row = linear_study.comparisonResults.Model == ...
    "Complex GMP DOMP parameter-matched";
if nnz(independent_row) ~= 1 || nnz(matched_row) ~= 1
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingParameterMatch', ...
        'The two parameter-matched linear rows must be present exactly once.');
end
target_active_params = ...
    linear_study.comparisonResults.NumRealParameters(independent_row);
if linear_study.comparisonResults.NumRealParameters(matched_row) ~= ...
        target_active_params
    error('run_fair_PNNN_vs_PNGMP_DOMP:ParameterMismatch', ...
        'Independent PN-IQ and parameter-matched GMP must have equal size.');
end

pnnn_study = runPNNNComparisonStudy( ...
    data.x, data.y, split, cfg, target_active_params);
study = combineComparisonStudies(linear_study, pnnn_study, split);
saveComparisonStudy(data, split, study, cfg, result_directory);
printComparisonSummary(study, split, result_directory);
diary off;
