% Script: run_domp_unit_test
% Verify the fixed DOMP selector against a direct pseudocode implementation
% for real, complex, mixed, and real multioutput deterministic fixtures.

clearvars;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'toolbox', 'domp'));

rng(2718, 'twister');
n_rows = 96;
n_columns = 14;
options = struct('columnTolerance', 1e-13, ...
    'correlationTolerance', 0, 'residualTolerance', 0, ...
    'lsTolerance', []);

X_real = randn(n_rows, n_columns);
y_real = X_real(:, [3, 8, 11]) * [1.2; -0.7; 0.35] + ...
    1e-5*randn(n_rows, 1);
checkFixture(X_real, y_real, 5, options);

X_complex = (randn(n_rows, n_columns) + ...
    1j*randn(n_rows, n_columns))/sqrt(2);
y_complex = X_complex(:, [2, 6, 13]) * ...
    [0.8-0.2j; -0.4+0.9j; 0.3+0.1j] + ...
    1e-5*(randn(n_rows, 1) + 1j*randn(n_rows, 1));
checkFixture(X_complex, y_complex, 5, options);

y_mixed = X_real(:, [4, 9]) * [0.9+0.3j; -0.5+0.7j] + ...
    1e-5*(randn(n_rows, 1) + 1j*randn(n_rows, 1));
checkFixture(X_real, y_mixed, 4, options);

y_multioutput = [real(y_mixed), imag(y_mixed)];
checkFixture(X_real, y_multioutput, 4, options);

fprintf('\nDOMP UNIT TEST: PASS\n');

function checkFixture(X, y, maximum_components, options)
[support, history] = selectDOMPSupport( ...
    X, y, maximum_components, options);
[support_repeat, history_repeat] = selectDOMPSupport( ...
    X, y, maximum_components, options);
reference_support = referenceDOMP( ...
    X, y, maximum_components, options);

assert(isequal(support, support_repeat));
assert(isequal(history.selectedIndex, history_repeat.selectedIndex));
assert(isequal(support, reference_support));
assert(numel(support) <= maximum_components);
assert(numel(unique(support)) == numel(support));
assert(all(support >= 1 & support <= size(X, 2)));
assert(all(isfinite(history.residualNorm)));
assert(all(isfinite(history.normalizedCorrelation)));
assert(all(isfinite(history.selectedColumnNorm)));
assert(all(isfinite(history.conditionNumberSelected)));
assert(all(isfinite(history.elapsedTime)));
assert(all(isfinite(history.maxOrthogonalityError)));
residual_increase = diff(history.residualNorm);
residual_scale = max(1, history.initialResidualNorm);
assert(all(residual_increase <= 1e-11*residual_scale));
assert(max(history.maxOrthogonalityError) <= 1e-10);
assert(history.selectedCount == numel(support));
assert(history.finalResidualNorm <= ...
    history.initialResidualNorm + 1e-11*residual_scale);
end

function support = referenceDOMP(X, y, maximum_components, options)
if isvector(y)
    y = y(:);
end
Z = X;
residual = y;
selected = false(size(X, 2), 1);
support = zeros(maximum_components, 1);
initial_norms = sqrt(sum(abs(X).^2, 1));
absolute_tolerance = options.columnTolerance * ...
    max(1, max(initial_norms));
n_selected = 0;

for iteration = 1:maximum_components
    norms = sqrt(sum(abs(Z).^2, 1));
    eligible = ~selected.' & norms > absolute_tolerance;
    if ~any(eligible)
        break;
    end
    Z(:, eligible) = Z(:, eligible) ./ norms(eligible);
    Z(:, ~eligible) = 0;
    eligible_indices = find(eligible);
    correlations = Z(:, eligible)' * residual;
    scores = sqrt(sum(abs(correlations).^2, 2));
    [best_score, local_index] = max(scores);
    if best_score <= options.correlationTolerance * ...
            max(1, norm(residual, 'fro'))
        break;
    end
    index = eligible_indices(local_index);
    q = Z(:, index);
    Z = Z - q * (q' * Z);
    Z(:, index) = 0;
    selected(index) = true;
    n_selected = n_selected + 1;
    support(n_selected) = index;
    active = support(1:n_selected);
    if isempty(options.lsTolerance)
        coefficients = lsqminnorm(X(:, active), y);
    else
        coefficients = lsqminnorm( ...
            X(:, active), y, options.lsTolerance);
    end
    residual = y - X(:, active)*coefficients;
    if norm(residual, 'fro') <= options.residualTolerance * ...
            max(1, norm(y, 'fro'))
        break;
    end
end
support = support(1:n_selected);
end
