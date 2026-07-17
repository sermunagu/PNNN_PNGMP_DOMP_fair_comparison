function results = main_sweep_and_comparison()
% main_sweep_and_comparison - Run the common sweep and select one budget.
% The signed sweep is resumed before prompting for a parameter count.
% The selected point is compared using only compatible sweep artifacts.

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
cfg = getFairDOMPComparisonConfig(projectRoot);
targets = cfg.sweep.parameterGrid;
targetText = mat2str(targets);
sweep = run_parameter_sweep(targets);
drawnow;

selectedParameters = NaN;
while ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ...
        ismember(selectedParameters, targets))
    userText = input([ ...
        'Review the sweep figures and enter the selected number of ' ...
        'parameters from ' targetText ': '], 's');
    selectedParameters = str2double(strtrim(userText));
    if ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ...
            ismember(selectedParameters, targets))
        fprintf('Invalid parameter budget. Choose a value from %s.\n', ...
            targetText);
    end
end

fprintf('Selected parameter budget: %d active real parameters.\n', ...
    selectedParameters);
results = run_selected_comparison(selectedParameters, sweep);
end
