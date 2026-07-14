function [yhat, yhat_phase_normalized] = predictIndependentIQGMP( ...
    x, rows, manager, support, fit, blockSize)
% predictIndependentIQGMP - Predict independent I/Q PN-GMP outputs by blocks.
% The same row phase is applied while constructing every block, and the
% saved structural feature map is reused without materializing full matrices.

if nargin < 6 || isempty(blockSize)
    blockSize = 8192;
end
rows = double(rows(:));
if isempty(rows) || any(~isfinite(rows)) || ...
        any(rows ~= floor(rows)) || any(rows < 1) || ...
        any(rows > numel(x))
    error('predictIndependentIQGMP:InvalidRows', ...
        'rows must contain valid signal indices.');
end
if ~isstruct(fit) || ~all(isfield(fit, ...
        {'coefficientsI','coefficientsQ','reduction'}))
    error('predictIndependentIQGMP:InvalidFit', ...
        'fit must contain independent I/Q coefficients and a reduction map.');
end
kept_indices = fit.reduction.keptIndices(:);
if numel(kept_indices) ~= fit.effectiveRealFeatures
    error('predictIndependentIQGMP:InvalidReduction', ...
        'The structural reduction map is inconsistent with the fit.');
end

blockSize = max(1, floor(double(blockSize)));
yhat_phase_normalized = complex(zeros(numel(rows), 1));
phase_rotation = complex(zeros(numel(rows), 1));
for first = 1:blockSize:numel(rows)
    last = min(first + blockSize - 1, numel(rows));
    block_rows = rows(first:last);
    [raw_features, details] = buildPhaseNormalizedIQRegressors( ...
        x, block_rows, manager, support);
    features = raw_features(:, kept_indices);
    prediction_I = features * fit.coefficientsI;
    prediction_Q = features * fit.coefficientsQ;
    yhat_phase_normalized(first:last) = ...
        prediction_I + 1j*prediction_Q;
    phase_rotation(first:last) = details.phaseRotation;
end
yhat = conj(phase_rotation) .* yhat_phase_normalized;
if any(~isfinite(yhat)) || any(~isfinite(yhat_phase_normalized))
    error('predictIndependentIQGMP:NonFinitePrediction', ...
        'The predictor produced NaN or Inf.');
end
end
