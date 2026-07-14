function [features, details] = buildPhaseNormalizedIQRegressors( ...
    x, rows, reg_manager, support_complex)
% buildPhaseNormalizedIQRegressors - Build explicit PN-IQ GMP features.
% Every canonical term uses the current-row phase for all delayed taps;
% noncanonical auxiliaries retain their exact normalized I/Q columns.

x = validateSignal(x);
rows = validateRows(rows, numel(x));
support_complex = validateSupport( ...
    support_complex, numel(reg_manager.regPopulation));
n_rows = numel(rows);
n_regressors = numel(support_complex);
phase_rotation = computePhaseNormGMPRotation(x, rows);
U_complex = buildGMPRegressorRows( ...
    x, rows, reg_manager, support_complex);
U_phase_normalized = phase_rotation .* U_complex;
regressors_I = zeros(n_rows, n_regressors);
regressors_Q = zeros(n_rows, n_regressors);
descriptors = repmat(factorizeGMPRegressor( ...
    reg_manager.regPopulation(support_complex(1)), support_complex(1)), ...
    n_regressors, 1);
canonical_I_errors = nan(n_regressors, 1);
canonical_Q_errors = nan(n_regressors, 1);
carrier_I0_errors = nan(n_regressors, 1);
carrier_Q0_maxima = nan(n_regressors, 1);

for local_index = 1:n_regressors
    population_index = support_complex(local_index);
    descriptor = factorizeGMPRegressor( ...
        reg_manager.regPopulation(population_index), population_index);
    descriptors(local_index) = descriptor;
    generic_I = real(U_phase_normalized(:, local_index));
    generic_Q = imag(U_phase_normalized(:, local_index));
    if descriptor.canonicalGMP
        carrier = x(wrapIndex( ...
            rows - descriptor.carrierLag, numel(x)));
        normalized_carrier = phase_rotation .* carrier;
        envelope_product = buildEnvelopeProduct( ...
            x, rows, descriptor, numel(x));
        regressors_I(:, local_index) = ...
            real(normalized_carrier) .* envelope_product;
        regressors_Q(:, local_index) = ...
            imag(normalized_carrier) .* envelope_product;
        if descriptor.QColumnStructurallyZero
            % The analytic identity r(n)x(n)=|x(n)| defines this column as
            % exactly zero; suppress roundoff from the generic complex path.
            regressors_Q(:, local_index) = 0;
        end
        canonical_I_errors(local_index) = relativeError( ...
            regressors_I(:, local_index), generic_I);
        if descriptor.QColumnStructurallyZero
            canonical_Q_errors(local_index) = norm(generic_Q) / ...
                max(1, norm(generic_I));
        else
            canonical_Q_errors(local_index) = relativeError( ...
                regressors_Q(:, local_index), generic_Q);
        end
        if descriptor.carrierLag == 0
            current_magnitude = abs(x(rows));
            carrier_I0_errors(local_index) = relativeError( ...
                real(normalized_carrier), current_magnitude);
            carrier_Q0_maxima(local_index) = ...
                max(abs(regressors_Q(:, local_index)));
        end
    else
        regressors_I(:, local_index) = generic_I;
        regressors_Q(:, local_index) = generic_Q;
    end
end

features = [regressors_I, regressors_Q];
SourceRegressorIndex = [support_complex; support_complex];
Component = [repmat("I", n_regressors, 1); ...
    repmat("Q", n_regressors, 1)];
Signature = strings(2*n_regressors, 1);
StructuralZero = false(2*n_regressors, 1);
CanonicalGMP = false(2*n_regressors, 1);
ExactAuxiliaryFallback = false(2*n_regressors, 1);
for local_index = 1:n_regressors
    descriptor = descriptors(local_index);
    Signature(local_index) = descriptor.iSignature;
    Signature(n_regressors + local_index) = descriptor.qSignature;
    StructuralZero(local_index) = descriptor.IColumnStructurallyZero;
    StructuralZero(n_regressors + local_index) = ...
        descriptor.QColumnStructurallyZero;
    CanonicalGMP([local_index, n_regressors + local_index]) = ...
        descriptor.canonicalGMP;
    ExactAuxiliaryFallback([local_index, n_regressors + local_index]) = ...
        ~descriptor.canonicalGMP;
end
RelationSign = ones(2*n_regressors, 1);
feature_metadata = table(SourceRegressorIndex, Component, Signature, ...
    StructuralZero, RelationSign, CanonicalGMP, ExactAuxiliaryFallback);

if any(~isfinite(features), 'all')
    error('buildPhaseNormalizedIQRegressors:NonFiniteOutput', ...
        'Explicit PN-IQ features contain NaN or Inf.');
end
details = struct();
details.phaseRotation = phase_rotation;
details.UPhaseNormalized = U_phase_normalized;
details.regressorsI = regressors_I;
details.regressorsQ = regressors_Q;
details.descriptors = descriptors;
details.featureMetadata = feature_metadata;
details.maxCanonicalIError = maxFinite(canonical_I_errors);
details.maxCanonicalQError = maxFinite(canonical_Q_errors);
details.maxCarrierI0Error = maxFinite(carrier_I0_errors);
details.maxCarrierQ0 = maxFinite(carrier_Q0_maxima);
details.canonicalCount = nnz([descriptors.canonicalGMP]);
details.noncanonicalCount = n_regressors - details.canonicalCount;
end

function envelope_product = buildEnvelopeProduct( ...
    x, rows, descriptor, n_signal)
envelope_product = ones(numel(rows), 1);
for term_index = 1:numel(descriptor.envelopeLags)
    magnitude = abs(x(wrapIndex( ...
        rows - descriptor.envelopeLags(term_index), n_signal)));
    envelope_product = envelope_product .* ...
        magnitude.^descriptor.envelopePowers(term_index);
end
end

function value = relativeError(actual, reference)
value = norm(actual - reference) / max(norm(reference), realmin);
end

function value = maxFinite(values)
values = values(isfinite(values));
if isempty(values)
    value = 0;
else
    value = max(values);
end
end

function idx = wrapIndex(idx, n_signal)
idx = mod(idx - 1, n_signal) + 1;
end

function x = validateSignal(x)
x = x(:);
if isempty(x) || ~isfloat(x) || any(~isfinite(x))
    error('buildPhaseNormalizedIQRegressors:InvalidSignal', ...
        'x must be a non-empty finite floating-point vector.');
end
end

function rows = validateRows(rows, upper_bound)
rows = double(rows(:));
if isempty(rows) || any(~isfinite(rows)) || ...
        any(rows ~= floor(rows)) || any(rows < 1) || ...
        any(rows > upper_bound)
    error('buildPhaseNormalizedIQRegressors:InvalidRows', ...
        'rows must contain valid signal indices.');
end
end

function support = validateSupport(support, upper_bound)
support = double(support(:));
if isempty(support) || any(~isfinite(support)) || ...
        any(support ~= floor(support)) || any(support < 1) || ...
        any(support > upper_bound) || ...
        numel(unique(support)) ~= numel(support)
    error('buildPhaseNormalizedIQRegressors:InvalidSupport', ...
        'support_complex must contain valid unique population indices.');
end
end
