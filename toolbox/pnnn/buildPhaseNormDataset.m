function [X_in, Y_out, r_vec] = buildPhaseNormDataset(x, y, M, orders, featMode)
% buildPhaseNormDataset - Build the phase-normalized PNNN dataset.
%
% This function creates input features, rotated two-channel targets, and
% phase rotations used by the offline PNNN training flow. Memory taps use
% periodic extension: x(n), x(n-1), ..., x(n-M).
%
% Inputs:
%   x, y - Modeled-block input and output signals under the local X/Y convention.
%   M, orders, featMode - Memory depth, nonlinear orders, and feature mode.
%
% Outputs:
%   X_in - D x N phase-normalized input feature matrix.
%   Y_out - 2 x N real/imag target matrix in the rotated frame.
%   r_vec - 1 x N phase-rotation vector for reconstruction.

    if nargin < 5 || isempty(featMode)
        featMode = "full";
    end

    x = x(:).';
    y = y(:).';

    if numel(x) ~= numel(y)
        error('x e y deben tener la misma longitud.');
    end
    if any(~isfinite(x)) || any(~isfinite(y))
        error('x o y contienen NaN/Inf.');
    end

    [X_in, r_vec] = buildPhaseNormInput(x, M, orders, featMode);
    N = numel(x);
    Y_out = zeros(2, N);

    for k = 1:N
        Y_out(:, k) = [real(r_vec(k)*y(k)); imag(r_vec(k)*y(k))];
    end
end
