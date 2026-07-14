function [support, history] = selectDOMPSupport( ...
    X, y, maxComponents, options)
% selectDOMPSupport - Select a sparse support using fixed DOMP.
% Candidate columns are Gram-Schmidt orthogonalized while the residual is
% recomputed by minimum-norm LS on the original design matrix each step.

if nargin < 4
    options = struct();
end
options = validateDOMPOptions(options);

if ~isnumeric(X) || ~isfloat(X) || ~ismatrix(X) || isempty(X) || ...
        any(~isfinite(X), 'all')
    error('selectDOMPSupport:InvalidDesign', ...
        'X must be a non-empty finite floating-point matrix.');
end
if ~isnumeric(y) || ~isfloat(y) || ~ismatrix(y) || isempty(y) || ...
        size(y, 1) ~= size(X, 1) || any(~isfinite(y), 'all')
    error('selectDOMPSupport:InvalidTarget', ...
        'y must be finite and have one row per observation in X.');
end
if isvector(y)
    y = y(:);
end
if ~isnumeric(maxComponents) || ~isreal(maxComponents) || ...
        ~isscalar(maxComponents) || ~isfinite(maxComponents) || ...
        maxComponents < 1 || maxComponents ~= floor(maxComponents)
    error('selectDOMPSupport:InvalidComponentCount', ...
        'maxComponents must be a positive integer scalar.');
end

n_candidates = size(X, 2);
maxComponents = min(double(maxComponents), n_candidates);
Z = X;
residual = y;
selected = false(n_candidates, 1);
support_buffer = zeros(maxComponents, 1);
selected_index = zeros(maxComponents, 1);
residual_norm = zeros(maxComponents, 1);
normalized_correlation = zeros(maxComponents, 1);
selected_column_norm = zeros(maxComponents, 1);
condition_number = zeros(maxComponents, 1);
elapsed_time = zeros(maxComponents, 1);
orthogonality_error = zeros(maxComponents, 1);

initial_column_norms = sqrt(sum(abs(X).^2, 1));
absolute_column_tolerance = options.columnTolerance * ...
    max(1, max(initial_column_norms));
initial_residual_norm = norm(y, 'fro');
absolute_residual_tolerance = options.residualTolerance * ...
    max(1, initial_residual_norm);
start_time = tic;
selected_count = 0;
stop_reason = "maximum_components";

for iteration = 1:maxComponents
    column_norms = sqrt(sum(abs(Z).^2, 1));
    eligible = ~selected.' & column_norms > absolute_column_tolerance;
    if ~any(eligible)
        stop_reason = "no_eligible_columns";
        break;
    end

    Z(:, eligible) = Z(:, eligible) ./ column_norms(eligible);
    Z(:, ~eligible) = 0;
    correlations = Z(:, eligible)' * residual;
    scores = sqrt(sum(abs(correlations).^2, 2));
    eligible_indices = find(eligible);
    [best_score, local_index] = max(scores);
    correlation_floor = options.correlationTolerance * ...
        max(1, norm(residual, 'fro'));
    if ~isfinite(best_score) || best_score <= correlation_floor
        stop_reason = "correlation_tolerance";
        break;
    end

    best_index = eligible_indices(local_index);
    q = Z(:, best_index);
    projections = q' * Z;
    Z = Z - q * projections;
    selected(best_index) = true;
    Z(:, best_index) = 0;

    selected_count = selected_count + 1;
    support_buffer(selected_count) = best_index;
    active_support = support_buffer(1:selected_count);
    X_selected = X(:, active_support);
    [coefficients, selected_condition] = solveSelectedLeastSquares( ...
        X_selected, y, options.lsTolerance);
    prediction = X_selected * coefficients;
    residual = y - prediction;

    remaining_norms = sqrt(sum(abs(Z).^2, 1));
    remaining = ~selected.' & ...
        remaining_norms > absolute_column_tolerance;
    if any(remaining)
        normalized_remaining = Z(:, remaining) ./ ...
            remaining_norms(remaining);
        orthogonality_error(selected_count) = ...
            max(abs(q' * normalized_remaining), [], 'all');
    end

    selected_index(selected_count) = best_index;
    residual_norm(selected_count) = norm(residual, 'fro');
    normalized_correlation(selected_count) = best_score;
    selected_column_norm(selected_count) = column_norms(best_index);
    condition_number(selected_count) = selected_condition;
    elapsed_time(selected_count) = toc(start_time);

    if residual_norm(selected_count) <= absolute_residual_tolerance
        stop_reason = "residual_tolerance";
        break;
    end
end

support = support_buffer(1:selected_count);
if isempty(support)
    error('selectDOMPSupport:EmptySupport', ...
        'DOMP did not find an eligible component.');
end
if numel(unique(support)) ~= numel(support)
    error('selectDOMPSupport:RepeatedIndex', ...
        'DOMP selected a repeated component.');
end

history = struct();
history.selectedIndex = selected_index(1:selected_count);
history.residualNorm = residual_norm(1:selected_count);
history.normalizedCorrelation = ...
    normalized_correlation(1:selected_count);
history.selectedColumnNorm = selected_column_norm(1:selected_count);
history.conditionNumberSelected = condition_number(1:selected_count);
history.elapsedTime = elapsed_time(1:selected_count);
history.maxOrthogonalityError = ...
    orthogonality_error(1:selected_count);
history.stopReason = stop_reason;
history.selectedCount = selected_count;
history.requestedComponents = maxComponents;
history.initialResidualNorm = initial_residual_norm;
history.finalResidualNorm = norm(residual, 'fro');
history.options = options;
history.targetColumns = size(y, 2);
end

function [coefficients, condition_number] = solveSelectedLeastSquares( ...
    X_selected, y, requested_tolerance)
% QR is the documented numerically equivalent LS solver for full-rank
% selected columns. Rank-deficient fixtures fall back to minimum-norm LS.
[Q, R] = qr(X_selected, 0);
if isempty(requested_tolerance)
    tolerance = max(size(X_selected))*eps(norm(R, 2));
else
    tolerance = requested_tolerance;
end
diagonal = abs(diag(R));
if isempty(diagonal) || any(diagonal <= tolerance)
    if isempty(requested_tolerance)
        coefficients = lsqminnorm(X_selected, y);
    else
        coefficients = lsqminnorm(X_selected, y, requested_tolerance);
    end
    condition_number = cond(X_selected);
else
    coefficients = R \ (Q' * y);
    condition_number = cond(R);
end
end
