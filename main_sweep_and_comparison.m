function results = main_sweep_and_comparison()
% main_sweep_and_comparison - Run the common sweep and select one budget.
% The signed sweep is resumed before prompting for a parameter count.
% The selected point is compared using only compatible sweep artifacts.

targets = 20:10:500;
sweep = run_parameter_sweep(targets);
drawnow;

selectedParameters = NaN;
while ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ...
        ismember(selectedParameters, targets))
    userText = input([ ...
        'Review the sweep figures and enter the selected number of ' ...
        'parameters (20:10:500): '], 's');
    selectedParameters = str2double(strtrim(userText));
    if ~(isscalar(selectedParameters) && isfinite(selectedParameters) && ...
            ismember(selectedParameters, targets))
        fprintf('Invalid parameter budget. Enter a value in 20:10:500.\n');
    end
end

fprintf('Selected parameter budget: %d active real parameters.\n', ...
    selectedParameters);
results = run_fair_PNNN_vs_PNGMP_DOMP(selectedParameters, sweep);
end
