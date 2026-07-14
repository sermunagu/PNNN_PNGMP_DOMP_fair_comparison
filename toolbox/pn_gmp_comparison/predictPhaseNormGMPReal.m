function [yhat, yhat_phase_normalized, details] = ...
    predictPhaseNormGMPReal( ...
    x, rows, rManagerGMP, fit, blockSize)
% predictPhaseNormGMPReal - Predict with the coupled real PN-GMP formulation.
% GMP regressors are evaluated from the original signal on each requested
% row, rotated together, predicted in real I/Q form, and rotated back to yhat.

x = x(:);
if isempty(x) || ~isnumeric(x) || any(~isfinite(x))
    error('predictPhaseNormGMPReal:InvalidInput', ...
        'x must be a non-empty finite numeric vector.');
end
if nargin < 2 || isempty(rows)
    rows = (1:numel(x)).';
end
if ~isnumeric(rows) || ~isreal(rows) || ~isvector(rows) || ...
        isempty(rows) || any(~isfinite(rows)) || ...
        any(rows ~= floor(rows)) || any(rows < 1) || any(rows > numel(x))
    error('predictPhaseNormGMPReal:InvalidRows', ...
        'rows must contain finite integer indices in [1, %d].', numel(x));
end
rows = double(rows(:));

if nargin < 4 || ~isstruct(fit) || ...
        ~isfield(fit, 'supportComplex') || ~isfield(fit, 'hReal')
    error('predictPhaseNormGMPReal:InvalidFit', ...
        'fit must contain supportComplex and hReal.');
end
supportComplex = fit.supportComplex(:);
hReal = fit.hReal(:);
if ~isreal(hReal) || any(~isfinite(hReal)) || ...
        numel(hReal) ~= 2*numel(supportComplex)
    error('predictPhaseNormGMPReal:InvalidCoefficients', ...
        'fit.hReal must contain two finite real values per complex regressor.');
end

if nargin < 5 || isempty(blockSize)
    blockSize = 8192;
end
if ~isscalar(blockSize) || ~isnumeric(blockSize) || ~isreal(blockSize) || ...
        ~isfinite(blockSize) || blockSize < 1 || blockSize ~= floor(blockSize)
    error('predictPhaseNormGMPReal:InvalidBlockSize', ...
        'blockSize must be a positive integer scalar.');
end
blockSize = min(double(blockSize), numel(rows));

nRows = numel(rows);
yhat = complex(zeros(nRows, 1));
yhat_phase_normalized = complex(zeros(nRows, 1));
phase_rotation_all = complex(zeros(nRows, 1));

for first = 1:blockSize:nRows
    last = min(first + blockSize - 1, nRows);
    local = first:last;
    blockRows = rows(local);

    U_block = buildGMPRegressorRows( ...
        x, blockRows, rManagerGMP, supportComplex);
    phase_rotation_block = computePhaseNormGMPRotation(x, blockRows);
    U_phase_normalized_block = phase_rotation_block .* U_block;
    [U_real_block, ~] = complexToCoupledReal( ...
        U_phase_normalized_block, []);
    yhat_real = U_real_block * hReal;

    nBlock = numel(blockRows);
    yhat_phase_normalized_block = complex( ...
        yhat_real(1:nBlock), yhat_real(nBlock+1:end));
    yhat_block = conj(phase_rotation_block) .* ...
        yhat_phase_normalized_block;

    yhat_phase_normalized(local) = yhat_phase_normalized_block;
    yhat(local) = yhat_block;
    phase_rotation_all(local) = phase_rotation_block;
end

if any(~isfinite(yhat)) || any(~isfinite(yhat_phase_normalized))
    error('predictPhaseNormGMPReal:NonFinitePrediction', ...
        'The PN-GMP prediction contains NaN or Inf values.');
end

details = struct();
details.rows = rows;
details.rotation = phase_rotation_all;
details.nActiveComplex = numel(supportComplex);
details.nActiveReal = numel(hReal);
details.primaryOutputField = 'yhat';
end
