function results = main_sweep_and_comparison(selectedParameters)
% main_sweep_and_comparison - Run the sweep and select one operating point.
% With no argument, the joint stabilization minimum-complexity criterion is
% used. An explicit signed-grid budget remains a supported manual override.

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
cfg = getFairDOMPComparisonConfig(projectRoot);

targets = cfg.sweep.parameterGrid;
sweep = run_parameter_sweep(targets);
drawnow;
automaticSelection = sweep.selection;
fprintf('Automatic selection: %s\n', automaticSelection.summarySentence);

if nargin < 1 || isempty(selectedParameters)
    selectedParameters = automaticSelection.selectedParameters;

else
    selectedParameters = double(selectedParameters);
    
    if ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ismember(selectedParameters, targets))
        error('main_sweep_and_comparison:InvalidManualOverride', 'Manual override must be one signed-grid budget from %s.', mat2str(targets));
    end
    
    fprintf(['Manual override requested: %d parameters; automatic ' 'selection remains %d parameters.\n'], selectedParameters, automaticSelection.selectedParameters);
end

fprintf('Selected parameter budget: %d active real parameters.\n', selectedParameters);

results = run_selected_comparison(selectedParameters, sweep);
results.selection = automaticSelection;
end