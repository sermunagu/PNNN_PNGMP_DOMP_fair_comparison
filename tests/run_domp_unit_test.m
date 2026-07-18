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
columnTolerance = 1e-13;

X_real = randn(n_rows, n_columns);
y_real = X_real(:, [3, 8, 11]) * [1.2; -0.7; 0.35] + ...
    1e-5*randn(n_rows, 1);
checkFixture(X_real, y_real, 5, columnTolerance);

X_complex = (randn(n_rows, n_columns) + ...
    1j*randn(n_rows, n_columns))/sqrt(2);
y_complex = X_complex(:, [2, 6, 13]) * ...
    [0.8-0.2j; -0.4+0.9j; 0.3+0.1j] + ...
    1e-5*(randn(n_rows, 1) + 1j*randn(n_rows, 1));
checkFixture(X_complex, y_complex, 5, columnTolerance);

y_mixed = X_real(:, [4, 9]) * [0.9+0.3j; -0.5+0.7j] + ...
    1e-5*(randn(n_rows, 1) + 1j*randn(n_rows, 1));
checkFixture(X_real, y_mixed, 4, columnTolerance);

y_multioutput = [real(y_mixed), imag(y_mixed)];
checkFixture(X_real, y_multioutput, 4, columnTolerance);

fprintf('\nDOMP UNIT TEST: PASS\n');

function checkFixture(X, y, maximum_components, columnTolerance)
support = selectDOMPSupport( ...
    X, y, maximum_components, columnTolerance);
support_repeat = selectDOMPSupport( ...
    X, y, maximum_components, columnTolerance);
reference_support = referenceDOMP( ...
    X, y, maximum_components, columnTolerance);

assert(isequal(support, support_repeat));
assert(isequal(support, reference_support));
assert(numel(support) == maximum_components);
assert(numel(unique(support)) == numel(support));
assert(all(support >= 1 & support <= size(X, 2)));
end

function support = referenceDOMP(X, y, maximum_components, columnTolerance)
if isvector(y)
    y = y(:);
end
Z = X;
residual = y;
selected = false(size(X, 2), 1);
support = zeros(maximum_components, 1);
initial_norms = sqrt(sum(abs(X).^2, 1));
absolute_tolerance = columnTolerance * ...
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
    [~, local_index] = max(scores);
    index = eligible_indices(local_index);
    q = Z(:, index);
    Z = Z - q * (q' * Z);
    Z(:, index) = 0;
    selected(index) = true;
    n_selected = n_selected + 1;
    support(n_selected) = index;
    active = support(1:n_selected);
    coefficients = lsqminnorm(X(:, active), y);
    residual = y - X(:, active)*coefficients;
end
support = support(1:n_selected);
end
