function stats = computePNNNNormalization(features, targets)
% computePNNNNormalization - Compute feature and target z-score statistics.
% Statistics are derived only from the rows supplied by the caller, allowing
% TRAIN-only selection and a separate FIT_POOL-only final refit.

if ~isnumeric(features) || ~isreal(features) || isempty(features) || ...
        any(~isfinite(features), 'all') || size(targets, 1) ~= size(features, 1) || ...
        size(targets, 2) ~= 2 || ~isreal(targets) || ...
        any(~isfinite(targets), 'all')
    error('computePNNNNormalization:InvalidData', ...
        'features and two-channel targets must be aligned finite real arrays.');
end

stats = struct();
stats.muX = mean(features, 1);
stats.sigmaX = std(features, 0, 1);
stats.sigmaX(stats.sigmaX == 0) = 1;
stats.muY = mean(targets, 1);
stats.sigmaY = std(targets, 0, 1);
stats.sigmaY(stats.sigmaY == 0) = 1;
end
