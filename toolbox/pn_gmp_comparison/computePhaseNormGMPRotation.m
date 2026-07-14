function phase_rotation = computePhaseNormGMPRotation(x, rows)
% computePhaseNormGMPRotation - Compute row-aligned PN-GMP phase rotations.
% Each selected row uses the phase of the modeled-block input at that row,
% with an exact unit rotation for zero-valued input samples.

x = x(:);
if isempty(x) || ~isnumeric(x) || ~isfloat(x) || any(~isfinite(x))
    error('computePhaseNormGMPRotation:InvalidInput', ...
        'x must be a non-empty finite floating-point vector.');
end

if nargin < 2 || isempty(rows)
    rows = (1:numel(x)).';
end
if ~isnumeric(rows) || ~isreal(rows) || ~isvector(rows) || ...
        isempty(rows) || any(~isfinite(rows)) || ...
        any(rows ~= floor(rows)) || any(rows < 1) || any(rows > numel(x))
    error('computePhaseNormGMPRotation:InvalidRows', ...
        'rows must contain finite integer indices in [1, %d].', numel(x));
end

rows = double(rows(:));
x_rows = x(rows);
phase_rotation = complex(ones(size(x_rows)));
nonzero = abs(x_rows) ~= 0;
phase_rotation(nonzero) = ...
    conj(x_rows(nonzero)) ./ abs(x_rows(nonzero));

if any(~isfinite(phase_rotation))
    error('computePhaseNormGMPRotation:NonFiniteRotation', ...
        'The computed phase rotation contains NaN or Inf values.');
end
end
