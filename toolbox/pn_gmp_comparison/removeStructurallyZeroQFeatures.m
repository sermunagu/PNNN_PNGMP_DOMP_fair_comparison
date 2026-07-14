function [features_reduced, reduction] = removeStructurallyZeroQFeatures( ...
    features, feature_metadata, tolerance)
% removeStructurallyZeroQFeatures - Remove only declared zero Q columns.
% No duplicate, opposite-sign, correlation, rank, or conditioning reduction
% is performed for the full and DOMP-reduced Independent PN-IQ models.

if nargin < 3 || isempty(tolerance)
    tolerance = 1e-12;
end
if ~isnumeric(features) || ~isreal(features) || isempty(features) || ...
        any(~isfinite(features), 'all')
    error('removeStructurallyZeroQFeatures:InvalidFeatures', ...
        'features must be a finite non-empty real matrix.');
end
required = {'Component','StructuralZero'};
if ~istable(feature_metadata) || height(feature_metadata) ~= size(features, 2) || ...
        ~all(ismember(required, feature_metadata.Properties.VariableNames))
    error('removeStructurallyZeroQFeatures:InvalidMetadata', ...
        'feature_metadata must identify every structural zero and component.');
end
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error('removeStructurallyZeroQFeatures:InvalidTolerance', ...
        'tolerance must be a positive finite scalar.');
end

zero_mask = logical(feature_metadata.StructuralZero(:));
if any(zero_mask & string(feature_metadata.Component) ~= "Q")
    error('removeStructurallyZeroQFeatures:UnexpectedZeroComponent', ...
        'Only structurally zero Q components may be removed.');
end
for index = find(zero_mask).'
    if norm(features(:, index), Inf) > tolerance
        error('removeStructurallyZeroQFeatures:DeclaredZeroMismatch', ...
            'Declared zero Q feature %d is numerically non-zero.', index);
    end
end

kept_indices = find(~zero_mask);
features_reduced = features(:, kept_indices);
raw_to_kept = zeros(size(features, 2), 1);
raw_to_kept(kept_indices) = (1:numel(kept_indices)).';
reduction = struct();
reduction.rawFeatureCount = size(features, 2);
reduction.effectiveFeatureCount = numel(kept_indices);
reduction.keptIndices = kept_indices;
reduction.zeroIndices = find(zero_mask);
reduction.duplicateIndices = zeros(0, 1);
reduction.oppositeIndices = zeros(0, 1);
reduction.rawToKept = raw_to_kept;
reduction.rawSign = ones(size(features, 2), 1);
reduction.structurallyZeroRemoved = nnz(zero_mask);
reduction.structuralDuplicatesRemoved = 0;
reduction.structuralOppositesRemoved = 0;
reduction.featureMetadata = feature_metadata;
reduction.tolerance = tolerance;
end
