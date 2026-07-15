function reduction = restrictIQFeatureReduction(baseReduction, localSupport)
% restrictIQFeatureReduction - Map selected real features to source metadata.
% The selected indices refer to the already reduced feature matrix and are
% stored together with their original feature indices for later prediction.

localSupport = double(localSupport(:));
reduction = baseReduction;
reduction.keptIndices = baseReduction.keptIndices(localSupport);
reduction.effectiveFeatureCount = numel(localSupport);
reduction.groupSupportFeatures = localSupport;
reduction.selectionMethod = 'DOMP';
end
