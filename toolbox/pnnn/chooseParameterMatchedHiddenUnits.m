function choice = chooseParameterMatchedHiddenUnits(inputDimension, targetParameters)
% chooseParameterMatchedHiddenUnits - Choose the closest one-layer PNNN size.
% The integer hidden width minimizes the absolute difference to the requested
% real-parameter target without claiming exact equality when none exists.

validateattributes(inputDimension, {'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(targetParameters, {'numeric'}, ...
    {'scalar','integer','positive','finite'});

continuous_width = (double(targetParameters) - 2) / ...
    (double(inputDimension) + 3);
candidates = unique(max(1, [floor(continuous_width), ...
    ceil(continuous_width)]));
parameter_counts = arrayfun(@(width) ...
    width*(double(inputDimension) + 3) + 2, candidates);
differences = abs(parameter_counts - double(targetParameters));
[~, best] = min(differences);

choice = struct();
choice.targetParameters = double(targetParameters);
choice.hiddenNeurons = double(candidates(best));
choice.realParameters = double(parameter_counts(best));
choice.absoluteDifference = choice.realParameters - ...
    choice.targetParameters;
choice.absoluteDifferenceMagnitude = abs(choice.absoluteDifference);
choice.relativeDifferencePercent = 100*choice.absoluteDifferenceMagnitude / ...
    choice.targetParameters;
end
