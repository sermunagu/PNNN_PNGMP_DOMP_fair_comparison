function [features, details] = buildUnnormalizedIQRegressors( ...
    x, rows, reg_manager, support_complex)
% buildUnnormalizedIQRegressors - Split an original complex GMP basis into I/Q.
% Only analytically real envelope-only columns are marked with structural Q
% zeros; no phase normalization or data-driven duplicate removal is applied.

x = x(:);
rows = double(rows(:));
support_complex = double(support_complex(:));
if isempty(x) || any(~isfinite(x)) || isempty(rows) || ...
        any(rows < 1) || any(rows > numel(x)) || ...
        any(rows ~= floor(rows)) || isempty(support_complex) || ...
        numel(unique(support_complex)) ~= numel(support_complex)
    error('buildUnnormalizedIQRegressors:InvalidInput', ...
        'Signal, rows, and support must be finite and valid.');
end

U_complex = buildGMPRegressorRows( ...
    x, rows, reg_manager, support_complex);
n_regressors = numel(support_complex);
features = [real(U_complex), imag(U_complex)];
SourceRegressorIndex = [support_complex; support_complex];
Component = [repmat("I", n_regressors, 1); ...
    repmat("Q", n_regressors, 1)];
Signature = strings(2*n_regressors, 1);
StructuralZero = false(2*n_regressors, 1);
RelationSign = ones(2*n_regressors, 1);
CanonicalGMP = false(2*n_regressors, 1);
ExactAuxiliaryFallback = false(2*n_regressors, 1);

for local_index = 1:n_regressors
    population_index = support_complex(local_index);
    regressor = reg_manager.regPopulation(population_index);
    descriptor = factorizeGMPRegressor(regressor, population_index);
    signature_root = "NO_PN:" + compose('%d', population_index);
    Signature(local_index) = signature_root + ":I";
    Signature(n_regressors + local_index) = signature_root + ":Q";
    CanonicalGMP([local_index, n_regressors + local_index]) = ...
        descriptor.canonicalGMP;
    ExactAuxiliaryFallback([local_index, ...
        n_regressors + local_index]) = ~descriptor.canonicalGMP;

    % A product containing envelopes only is analytically real before PN.
    is_envelope_only = isempty(regressor.X) && ...
        isempty(regressor.Xconj) && ~isempty(regressor.Xenv);
    StructuralZero(n_regressors + local_index) = is_envelope_only;

    is_current_linear = descriptor.canonicalGMP && ...
        descriptor.carrierLag == 0 && isempty(descriptor.envelopeLags);
    is_current_conjugate = descriptor.auxiliaryType == "conjugate";
    if is_current_linear || is_current_conjugate
        Signature(local_index) = "NO_PN:CURRENT_LINEAR:I";
        Signature(n_regressors + local_index) = ...
            "NO_PN:CURRENT_LINEAR:Q";
        if is_current_conjugate
            RelationSign(n_regressors + local_index) = -1;
        end
    end
end

feature_metadata = table(SourceRegressorIndex, Component, Signature, ...
    StructuralZero, RelationSign, CanonicalGMP, ExactAuxiliaryFallback);
if any(~isfinite(features), 'all')
    error('buildUnnormalizedIQRegressors:NonFiniteOutput', ...
        'The unnormalized I/Q basis contains NaN or Inf.');
end
details = struct();
details.UComplex = U_complex;
details.featureMetadata = feature_metadata;
details.structurallyZeroQCount = nnz(StructuralZero);
details.rawFeatureCount = size(features, 2);
end
