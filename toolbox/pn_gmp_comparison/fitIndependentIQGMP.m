function fit = fitIndependentIQGMP( ...
    features, y_phase_normalized, lambda, reduction, model_name)
% fitIndependentIQGMP - Fit two real outputs on one normalized feature basis.
% LS uses a common minimum-norm rank tolerance; ridge applies the same
% lambda independently to phase-normalized I and Q targets.

if nargin < 3 || isempty(lambda)
    lambda = 0;
end
if nargin < 4
    reduction = struct();
end
if nargin < 5 || isempty(model_name)
    model_name = "Independent PN-IQ-GMP";
end
if ~isnumeric(features) || ~isreal(features) || isempty(features) || ...
        any(~isfinite(features), 'all')
    error('fitIndependentIQGMP:InvalidFeatures', ...
        'features must be a finite non-empty real matrix.');
end
y_phase_normalized = y_phase_normalized(:);
if numel(y_phase_normalized) ~= size(features, 1) || ...
        any(~isfinite(y_phase_normalized))
    error('fitIndependentIQGMP:InvalidTarget', ...
        'y_phase_normalized must contain one finite value per feature row.');
end
if ~isscalar(lambda) || ~isreal(lambda) || ...
        ~isfinite(lambda) || lambda < 0
    error('fitIndependentIQGMP:InvalidLambda', ...
        'lambda must be a finite non-negative scalar.');
end

feature_norms = sqrt(sum(features.^2, 1)).';
if any(feature_norms == 0)
    error('fitIndependentIQGMP:ZeroFeature', ...
        'Structural reduction must remove zero features before fitting.');
end
features_normalized = features ./ feature_norms.';
gram_normalized = features_normalized.' * features_normalized;
target_I = real(y_phase_normalized);
target_Q = imag(y_phase_normalized);
rhs_I = features_normalized.' * target_I;
rhs_Q = features_normalized.' * target_Q;
spectral_norm = norm(features_normalized, 2);
rank_tolerance = max(size(features_normalized)) * eps(spectral_norm);
effective_rank = rank(features_normalized, rank_tolerance);

if lambda == 0
    coefficients_I_normalized = lsqminnorm( ...
        features_normalized, target_I, rank_tolerance);
    coefficients_Q_normalized = lsqminnorm( ...
        features_normalized, target_Q, rank_tolerance);
    solver_name = "LS";
else
    regularized_gram = gram_normalized + ...
        lambda*eye(size(gram_normalized));
    coefficients_I_normalized = regularized_gram \ rhs_I;
    coefficients_Q_normalized = regularized_gram \ rhs_Q;
    if lambda == 1e-3
        solver_name = "ridge 1e-3";
    elseif lambda == 1e-4
        solver_name = "ridge 1e-4";
    else
        solver_name = "ridge " + compose('%g', lambda);
    end
end
coefficients_I = coefficients_I_normalized ./ feature_norms;
coefficients_Q = coefficients_Q_normalized ./ feature_norms;
if any(~isfinite([coefficients_I; coefficients_Q]))
    error('fitIndependentIQGMP:NonFiniteCoefficients', ...
        'Fitted coefficients contain NaN or Inf.');
end

fit = struct();
fit.modelName = char(model_name);
fit.lambda = double(lambda);
fit.solver = char(solver_name);
fit.coefficientsI = coefficients_I;
fit.coefficientsQ = coefficients_Q;
fit.coefficientsINormalized = coefficients_I_normalized;
fit.coefficientsQNormalized = coefficients_Q_normalized;
fit.featureNorms = feature_norms;
fit.gramNormalized = gram_normalized;
fit.gramConditionNumber = cond(gram_normalized);
fit.gramConditionNumberI = fit.gramConditionNumber;
fit.gramConditionNumberQ = fit.gramConditionNumber;
fit.rankTolerance = rank_tolerance;
fit.effectiveRank = effective_rank;
fit.effectiveRealFeatures = size(features, 2);
fit.realParameters = 2*size(features, 2);
fit.coefficientMemoryBytes = 8*fit.realParameters;
fit.reduction = reduction;
end
