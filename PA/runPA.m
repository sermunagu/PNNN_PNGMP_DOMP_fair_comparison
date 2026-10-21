function results = runPA(selectedParameters)
% runPA - Run the PA sweep and compare one operating point.
% With no argument, the joint stabilization minimum-complexity criterion is
% used. An explicit signed-grid budget remains a supported manual override.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
cfg = getComparisonConfig(projectRoot);

targets = cfg.sweep.parameterGrid;
sweep = runPASweep(targets);
drawnow;
automaticSelection = sweep.selection;
fprintf('Automatic selection: %s\n', automaticSelection.summarySentence);

if nargin < 1 || isempty(selectedParameters)
    selectedParameters = automaticSelection.selectedParameters;

else
    selectedParameters = double(selectedParameters);
    
    if ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ismember(selectedParameters, targets))
        error('runPA:InvalidManualOverride', 'Manual override must be one signed-grid budget from %s.', mat2str(targets));
    end
    
    fprintf(['Manual override requested: %d parameters; automatic ' 'selection remains %d parameters.\n'], selectedParameters, automaticSelection.selectedParameters);
end

fprintf('Selected parameter budget: %d active real parameters.\n', selectedParameters);

results = comparePA(selectedParameters, sweep);
results.selection = automaticSelection;
end
