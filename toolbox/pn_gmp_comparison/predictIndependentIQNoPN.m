function yhat = predictIndependentIQNoPN( ...
    x, rows, manager, support, fit, block_size)
% predictIndependentIQNoPN - Predict the independent I/Q control without PN.
% The frozen complex DOMP support and structural reduction map are applied
% by blocks, with no phase rotation or phase-restoration operation.

if nargin < 6 || isempty(block_size)
    block_size = 8192;
end
rows = double(rows(:));
if isempty(rows) || any(~isfinite(rows)) || ...
        any(rows ~= floor(rows)) || any(rows < 1) || any(rows > numel(x))
    error('predictIndependentIQNoPN:InvalidRows', ...
        'rows must contain valid signal indices.');
end
required = {'coefficientsI','coefficientsQ','reduction', ...
    'effectiveRealFeatures'};
if ~isstruct(fit) || ~all(isfield(fit, required))
    error('predictIndependentIQNoPN:InvalidFit', ...
        'fit does not satisfy the independent I/Q contract.');
end

yhat = complex(zeros(numel(rows), 1));
block_size = max(1, floor(double(block_size)));
for first = 1:block_size:numel(rows)
    last = min(first + block_size - 1, numel(rows));
    block_rows = rows(first:last);
    raw_features = buildUnnormalizedIQRegressors( ...
        x, block_rows, manager, support);
    features = raw_features(:, fit.reduction.keptIndices);
    yhat(first:last) = features*fit.coefficientsI + ...
        1j*(features*fit.coefficientsQ);
end
if any(~isfinite(yhat))
    error('predictIndependentIQNoPN:NonFinitePrediction', ...
        'The predictor produced NaN or Inf.');
end
end
