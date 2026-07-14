function [features_reduced, reduction] = reduceStructuralFeatures( ...
    features, feature_metadata, tolerance)
% reduceStructuralFeatures - Remove only analytically declared zeros and copies.
% Signatures define structural equivalence; numerical checks validate those
% declarations but never discover reductions from conditioning or correlation.

if nargin < 3 || isempty(tolerance)
    tolerance = 1e-12;
end
if ~isnumeric(features) || ~isreal(features) || isempty(features) || ...
        any(~isfinite(features), 'all')
    error('reduceStructuralFeatures:InvalidFeatures', ...
        'features must be a finite non-empty real matrix.');
end
required_variables = ["Signature", "StructuralZero", "RelationSign"];
if ~istable(feature_metadata) || height(feature_metadata) ~= size(features, 2) || ...
        ~all(ismember(required_variables, ...
        string(feature_metadata.Properties.VariableNames)))
    error('reduceStructuralFeatures:InvalidMetadata', ...
        'feature_metadata must describe every raw feature.');
end
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error('reduceStructuralFeatures:InvalidTolerance', ...
        'tolerance must be a positive finite scalar.');
end

n_raw = size(features, 2);
zero_mask = logical(feature_metadata.StructuralZero(:));
for index = find(zero_mask).'
    absolute_tolerance = tolerance * max(1, norm(features(:, index), Inf));
    if norm(features(:, index), Inf) > absolute_tolerance
        error('reduceStructuralFeatures:DeclaredZeroMismatch', ...
            'Structurally zero feature %d is numerically nonzero.', index);
    end
end

kept_indices = zeros(0, 1);
duplicate_indices = zeros(0, 1);
opposite_indices = zeros(0, 1);
raw_to_kept = zeros(n_raw, 1);
raw_sign = ones(n_raw, 1);
kept_declared_sign = zeros(0, 1);
signature_to_position = containers.Map('KeyType', 'char', ...
    'ValueType', 'double');
for index = 1:n_raw
    if zero_mask(index)
        continue;
    end
    signature = char(feature_metadata.Signature(index));
    declared_sign = double(feature_metadata.RelationSign(index));
    if ~any(declared_sign == [-1, 1])
        error('reduceStructuralFeatures:InvalidRelationSign', ...
            'RelationSign must be +1 or -1.');
    end
    if ~isKey(signature_to_position, signature)
        kept_indices(end+1, 1) = index; %#ok<AGROW>
        position = numel(kept_indices);
        signature_to_position(signature) = position;
        raw_to_kept(index) = position;
        raw_sign(index) = 1;
        kept_declared_sign(position, 1) = declared_sign;
        continue;
    end
    position = signature_to_position(signature);
    reference_index = kept_indices(position);
    relative_sign = declared_sign / kept_declared_sign(position);
    expected = relative_sign * features(:, reference_index);
    relative_error = norm(features(:, index) - expected) / ...
        max(1, norm(expected));
    if relative_error > tolerance
        error('reduceStructuralFeatures:SignatureMismatch', ...
            ['Features %d and %d share a structural signature but differ ' ...
             'numerically (relative error %.3e).'], ...
            reference_index, index, relative_error);
    end
    raw_to_kept(index) = position;
    raw_sign(index) = relative_sign;
    if relative_sign == 1
        duplicate_indices(end+1, 1) = index; %#ok<AGROW>
    else
        opposite_indices(end+1, 1) = index; %#ok<AGROW>
    end
end

features_reduced = features(:, kept_indices);
reduction = struct();
reduction.rawFeatureCount = n_raw;
reduction.effectiveFeatureCount = numel(kept_indices);
reduction.keptIndices = kept_indices;
reduction.zeroIndices = find(zero_mask);
reduction.duplicateIndices = duplicate_indices;
reduction.oppositeIndices = opposite_indices;
reduction.rawToKept = raw_to_kept;
reduction.rawSign = raw_sign;
reduction.keptDeclaredSign = kept_declared_sign;
reduction.structurallyZeroRemoved = nnz(zero_mask);
reduction.structuralDuplicatesRemoved = numel(duplicate_indices);
reduction.structuralOppositesRemoved = numel(opposite_indices);
reduction.featureMetadata = feature_metadata;
reduction.tolerance = tolerance;
end
