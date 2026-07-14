function [summary, analysis] = analyzeRegressorStructure( ...
    population, support)
% analyzeRegressorStructure - Summarize the retained GMP and PN-IQ structure.
% Regressor metadata is derived from X, Xconj, and Xenv, including both
% noncanonical auxiliary terms and exact I/Q structural identities.

n_population = numel(population);
support = double(support(:));
if n_population == 0 || isempty(support) || ...
        any(~isfinite(support)) || any(support ~= floor(support)) || ...
        any(support < 1) || any(support > n_population) || ...
        numel(unique(support)) ~= numel(support)
    error('analyzeRegressorStructure:InvalidInputs', ...
        'Population and support must contain valid unique regressors.');
end

descriptors = repmat(factorizeGMPRegressor(population(1), 1), ...
    n_population, 1);
for index = 1:n_population
    descriptors(index) = factorizeGMPRegressor(population(index), index);
end

population_mask = true(n_population, 1);
support_mask = ismember((1:n_population).', support);
population_counts = summarize(descriptors, population_mask);
support_counts = summarize(descriptors, support_mask);

Scope = ["Population"; "Frozen complex support"];
PopulationCount = repmat(n_population, 2, 1);
SupportCount = repmat(numel(support), 2, 1);
CanonicalCount = [population_counts.canonical; support_counts.canonical];
NoncanonicalCount = [population_counts.noncanonical; ...
    support_counts.noncanonical];
StructurallyRealCount = [population_counts.structurallyReal; ...
    support_counts.structurallyReal];
StructurallyZeroQCount = [population_counts.zeroQ; support_counts.zeroQ];
UniqueIQFeatureCount = [population_counts.uniqueIQ; support_counts.uniqueIQ];
summary = table(Scope, PopulationCount, SupportCount, CanonicalCount, ...
    NoncanonicalCount, StructurallyRealCount, StructurallyZeroQCount, ...
    UniqueIQFeatureCount);

regressor_indices = (1:n_population).';
canonical_mask = [descriptors.canonicalGMP].';
analysis = struct();
analysis.descriptors = descriptors;
analysis.population = population_counts;
analysis.support = support_counts;
analysis.supportComplex = support;
analysis.auxiliaryPopulationIndices = regressor_indices(~canonical_mask);
analysis.auxiliarySupportIndices = regressor_indices( ...
    ~canonical_mask & support_mask);
end

function counts = summarize(descriptors, mask)
selected = find(mask(:));
canonical = false(numel(selected), 1);
structurally_real = false(numel(selected), 1);
zero_q = false(numel(selected), 1);
iq_signatures = strings(0, 1);

for local_index = 1:numel(selected)
    descriptor = descriptors(selected(local_index));
    canonical(local_index) = descriptor.canonicalGMP;
    structurally_real(local_index) = ...
        descriptor.structurallyRealAfterPhaseNormalization;
    zero_q(local_index) = descriptor.QColumnStructurallyZero;
    iq_signatures(end+1, 1) = descriptor.iSignature; %#ok<AGROW>
    if ~descriptor.QColumnStructurallyZero
        iq_signatures(end+1, 1) = descriptor.qSignature; %#ok<AGROW>
    end
end

counts = struct();
counts.count = numel(selected);
counts.canonical = nnz(canonical);
counts.noncanonical = numel(selected) - counts.canonical;
counts.structurallyReal = nnz(structurally_real);
counts.zeroQ = nnz(zero_q);
counts.uniqueIQ = numel(unique(iq_signatures, 'stable'));
counts.rawIQ = 2*numel(selected);
end
