function fit = fitComplexGMPGrid( ...
    U_population, target, support, lambda_grid, normU_full)
% fitComplexGMPGrid - Fit one complex GMP support over a lambda grid.
% Every candidate uses the same identification-domain column normalization;
% lambda zero uses minimum-norm LS and positive lambdas use ridge equations.

if nargin < 5
    normU_full = [];
end
if ~isnumeric(U_population) || isempty(U_population) || ...
        any(~isfinite(U_population), 'all')
    error('fitComplexGMPGrid:InvalidDesign', ...
        'U_population must be a finite non-empty matrix.');
end
target = target(:);
support = double(support(:));
lambda_grid = double(lambda_grid(:));
if numel(target) ~= size(U_population, 1) || any(~isfinite(target)) || ...
        isempty(support) || any(support < 1) || ...
        any(support > size(U_population, 2)) || ...
        numel(unique(support)) ~= numel(support) || ...
        isempty(lambda_grid) || any(~isfinite(lambda_grid)) || ...
        any(lambda_grid < 0)
    error('fitComplexGMPGrid:InvalidInput', ...
        'Target, support, or lambda grid is invalid.');
end

if isempty(normU_full)
    normU_full = sqrt(sum(abs(U_population).^2, 1)).';
else
    normU_full = double(normU_full(:));
end
if numel(normU_full) ~= size(U_population, 2) || ...
        any(~isfinite(normU_full)) || any(normU_full <= 0)
    error('fitComplexGMPGrid:InvalidNorms', ...
        'normU_full must contain one positive norm per population column.');
end

selected_norms = normU_full(support);
U_normalized = U_population(:, support) ./ selected_norms.';
gram = U_normalized' * U_normalized;
rhs = U_normalized' * target;
rank_tolerance = max(size(U_normalized))*eps(norm(U_normalized, 2));
n_lambdas = numel(lambda_grid);
normalized_coefficients = complex(zeros(numel(support), n_lambdas));
coefficients = complex(zeros(numel(support), n_lambdas));
for index = 1:n_lambdas
    lambda = lambda_grid(index);
    if lambda == 0
        normalized_coefficients(:, index) = lsqminnorm( ...
            U_normalized, target, rank_tolerance);
    else
        normalized_coefficients(:, index) = ...
            (gram + lambda*eye(numel(support))) \ rhs;
    end
    coefficients(:, index) = ...
        normalized_coefficients(:, index) ./ selected_norms;
end

fit = struct();
fit.support = support;
fit.lambdaGrid = lambda_grid;
fit.coefficients = coefficients;
fit.normalizedCoefficients = normalized_coefficients;
fit.normU = normU_full;
fit.gramNormalized = gram;
fit.conditionNumber = cond(gram);
fit.rankTolerance = rank_tolerance;
fit.selectionMethod = 'DOMP';
end
