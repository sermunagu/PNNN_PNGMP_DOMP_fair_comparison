function metrics = evaluateIndependentIQFits( ...
    x, y, rows, reg_manager, support_complex, fits, ...
    representation, block_size)
% evaluateIndependentIQFits - Evaluate several independent-I/Q fits per block.
% One feature block is shared by all solvers, avoiding repeated basis builds;
% phase restoration is applied only for the explicit PN-IQ representation.

if nargin < 8 || isempty(block_size)
    block_size = 8192;
end
representation = lower(string(representation));
if ~any(representation == ["phase_normalized", "no_phase_normalization"])
    error('evaluateIndependentIQFits:InvalidRepresentation', ...
        'Unknown independent-I/Q representation.');
end
if ~iscell(fits) || isempty(fits)
    error('evaluateIndependentIQFits:InvalidFits', ...
        'fits must be a non-empty cell array.');
end
x = x(:);
y = y(:);
rows = double(rows(:));
n_solvers = numel(fits);
coefficients_I = zeros(fits{1}.effectiveRealFeatures, n_solvers);
coefficients_Q = zeros(fits{1}.effectiveRealFeatures, n_solvers);
for solver_index = 1:n_solvers
    if fits{solver_index}.effectiveRealFeatures ~= size(coefficients_I, 1) || ...
            ~isequal(fits{solver_index}.reduction.keptIndices, ...
            fits{1}.reduction.keptIndices)
        error('evaluateIndependentIQFits:InconsistentFits', ...
            'All solvers must share the same feature basis and reduction.');
    end
    coefficients_I(:, solver_index) = fits{solver_index}.coefficientsI;
    coefficients_Q(:, solver_index) = fits{solver_index}.coefficientsQ;
end

predictions_all = complex(zeros(numel(rows), n_solvers));
for first = 1:block_size:numel(rows)
    last = min(first + block_size - 1, numel(rows));
    block_rows = rows(first:last);
    if representation == "phase_normalized"
        [raw_features, details] = buildPhaseNormalizedIQRegressors( ...
            x, block_rows, reg_manager, support_complex);
        phase_rotation = details.phaseRotation;
    else
        [raw_features, ~] = buildUnnormalizedIQRegressors( ...
            x, block_rows, reg_manager, support_complex);
        phase_rotation = ones(numel(block_rows), 1);
    end
    features = raw_features(:, fits{1}.reduction.keptIndices);
    prediction_phase = features*coefficients_I + ...
        1j*(features*coefficients_Q);
    predictions = conj(phase_rotation) .* prediction_phase;
    predictions_all(first:last, :) = predictions;
end

metrics = struct();
target = y(rows);
metrics.NMSE_dB = zeros(1, n_solvers);
for solver_index = 1:n_solvers
    metrics.NMSE_dB(solver_index) = nmseComplexDb( ...
        target, predictions_all(:, solver_index));
end
metrics.errorEnergy = sum(abs(predictions_all - target).^2, 1);
metrics.targetEnergy = sum(abs(target).^2);
metrics.rows = numel(rows);
metrics.representation = char(representation);
end
