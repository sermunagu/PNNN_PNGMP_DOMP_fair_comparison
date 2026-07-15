function [comparison, details] = countModelOperations( ...
    population, support, iqReduction)
% countModelOperations - Count analytical real operations for retained models.
% Stage counts use explicit literal or reuse-aware schedules; they are not
% weighted FLOPs, execution-time estimates, or hardware measurements.

descriptors = describePopulation(population);
support = validateSupport(support, numel(population));
support_descriptors = descriptors(support);
S = numel(support);
K = nnz([support_descriptors.structurallyRealAfterPhaseNormalization]);

standard_dependencies = buildStandardDependencies(support_descriptors);
standard_generation = standardGenerationCost( ...
    support_descriptors, standard_dependencies);
structure_dependencies = buildStructureDependencies( ...
    support_descriptors, standard_dependencies);
structure_generation = structureGenerationCost( ...
    support_descriptors, structure_dependencies);

iq_signatures = collectIQSignatures(support_descriptors);
F_iq = numel(unique(iq_signatures, 'stable'));
if nargin >= 3 && ~isempty(iqReduction)
    F_iq = iqReduction.effectiveFeatureCount;
end

standard_M = standard_generation.M;
standard_A = standard_generation.A;
structure_M = structure_generation.M;
structure_A = structure_generation.A;
n_rotated_carriers = numel(structure_dependencies.nonzeroCarrierLags);
has_aux_conjugate = any(arrayfun(@(item) ...
    item.auxiliaryType == "conjugate", support_descriptors));
n_rotated_exact_auxiliaries = double(has_aux_conjugate) + ...
    nnz(arrayfun(@(item) item.family == ...
    "unsupported/noncanonical term", support_descriptors));
rotation_entities = n_rotated_carriers + n_rotated_exact_auxiliaries;

Model = [ ...
    "Complex GMP"; ...
    "Coupled PN-IQ-GMP materialized"; ...
    "Coupled PN-IQ-GMP structure-aware"; ...
    "Independent PN-IQ-GMP"];
Schedule = [ ...
    "complex literal"; ...
    "materialized coupled"; ...
    "reuse-aware analytical"; ...
    "reuse-aware analytical"];
n_models = numel(Model);

FeatureGenerationMultiplications = zeros(n_models, 1);
FeatureGenerationAdditions = zeros(n_models, 1);
CoefficientMultiplications = zeros(n_models, 1);
CoefficientProductAdditions = zeros(n_models, 1);
AccumulationAdditions = zeros(n_models, 1);
PhaseNormalizationMultiplications = zeros(n_models, 1);
PhaseNormalizationAdditions = zeros(n_models, 1);
PhaseRestorationMultiplications = zeros(n_models, 1);
PhaseRestorationAdditions = zeros(n_models, 1);
Divisions = zeros(n_models, 1);
SquareRoots = zeros(n_models, 1);
Comparisons = zeros(n_models, 1);
EffectiveRealFeatures = [2*S; 2*S; 2*S; F_iq];
RealParameters = [2*S; 2*S; 2*S; 2*F_iq];

% Complex GMP: literal complex products and complex accumulation.
FeatureGenerationMultiplications(1) = standard_M;
FeatureGenerationAdditions(1) = standard_A;
CoefficientMultiplications(1) = 4*S;
CoefficientProductAdditions(1) = 2*S;
AccumulationAdditions(1) = 2*max(S - 1, 0);
SquareRoots(1) = standard_dependencies.nEnvelopeLags;

% Materialized coupled form, including explicit row rotations.
FeatureGenerationMultiplications(2) = standard_M;
FeatureGenerationAdditions(2) = standard_A;
CoefficientMultiplications(2) = 4*S;
AccumulationAdditions(2) = 2*max(2*S - 1, 0);
PhaseNormalizationMultiplications(2) = 4 + 4*S;
PhaseNormalizationAdditions(2) = 2 + 2*S;
PhaseRestorationMultiplications(2) = 4;
PhaseRestorationAdditions(2) = 2;
Divisions(2) = 2;
SquareRoots(2) = standard_dependencies.nEnvelopeLags + 2;
Comparisons(2) = 1;

% Structure-aware coupled and independent I/Q share arithmetic, not memory.
for row = [3, 4]
    FeatureGenerationMultiplications(row) = structure_M;
    FeatureGenerationAdditions(row) = structure_A;
    PhaseNormalizationMultiplications(row) = 4*rotation_entities;
    PhaseNormalizationAdditions(row) = 2*rotation_entities;
    PhaseRestorationMultiplications(row) = 4;
    PhaseRestorationAdditions(row) = 2;
    Divisions(row) = 2;
    SquareRoots(row) = structure_dependencies.nEnvelopeLags;
    Comparisons(row) = 1;
end
CoefficientMultiplications(3) = 4*(S - K) + 2*K;
CoefficientProductAdditions(3) = 2*(S - K);
AccumulationAdditions(3) = 2*max(S - 1, 0);
CoefficientMultiplications(4) = 2*F_iq;
AccumulationAdditions(4) = 2*max(F_iq - 1, 0);

CoefficientMemoryBytes = 8*RealParameters;
TotalRealMultiplicationsPerSample = ...
    FeatureGenerationMultiplications + CoefficientMultiplications + ...
    PhaseNormalizationMultiplications + PhaseRestorationMultiplications;
TotalRealAdditionsPerSample = FeatureGenerationAdditions + ...
    CoefficientProductAdditions + AccumulationAdditions + ...
    PhaseNormalizationAdditions + PhaseRestorationAdditions;

comparison = table(Model, Schedule, EffectiveRealFeatures, ...
    RealParameters, CoefficientMemoryBytes, ...
    FeatureGenerationMultiplications, FeatureGenerationAdditions, ...
    CoefficientMultiplications, CoefficientProductAdditions, ...
    AccumulationAdditions, PhaseNormalizationMultiplications, ...
    PhaseNormalizationAdditions, PhaseRestorationMultiplications, ...
    PhaseRestorationAdditions, Divisions, SquareRoots, Comparisons, ...
    TotalRealMultiplicationsPerSample, TotalRealAdditionsPerSample);

numeric_counts = comparison{:, 3:end};
if any(~isfinite(numeric_counts), 'all') || ...
        any(numeric_counts < 0, 'all') || ...
        any(numeric_counts ~= floor(numeric_counts), 'all')
    error('countModelOperations:InvalidCount', ...
        'All operation and memory counts must be non-negative integers.');
end
if S == 100 && K == 17 && F_iq == 172 && ( ...
        comparison.TotalRealMultiplicationsPerSample(1) ~= 779 || ...
        comparison.TotalRealAdditionsPerSample(1) ~= 413 || ...
        comparison.TotalRealMultiplicationsPerSample(2) ~= 1187 || ...
        comparison.TotalRealAdditionsPerSample(2) ~= 617 || ...
        comparison.TotalRealMultiplicationsPerSample(3) ~= 774 || ...
        comparison.TotalRealAdditionsPerSample(3) ~= 403 || ...
        comparison.TotalRealMultiplicationsPerSample(4) ~= 774 || ...
        comparison.TotalRealAdditionsPerSample(4) ~= 403)
    error('countModelOperations:BaselineRegression', ...
        'Previously validated operation totals changed unexpectedly.');
end

details = struct();
details.supportComplex = support;
details.activeComplexRegressors = S;
details.phaseNormalizedStructurallyRealCount = K;
details.phaseNormalizedStructurallyRealIndices = support( ...
    [support_descriptors.structurallyRealAfterPhaseNormalization].');
details.standardDependencies = standard_dependencies;
details.structureAwareDependencies = structure_dependencies;
details.standardGenerationCost = standard_generation;
details.structureAwareGenerationCost = structure_generation;
details.independentIQEffectiveFeatures = F_iq;
details.assumptions = struct( ...
    'generalComplexMultiply', '4M+2A', ...
    'complexByRealMultiply', '2M', ...
    'complexAccumulation', '2A', ...
    'complexModulus', '2M+1A+1SQRT', ...
    'phaseDivision', '2DIV', ...
    'coefficientScalarBytes', 8, ...
    'weightedFlops', false, ...
    'hardwareSpeedupDemonstrated', false);
end

function descriptors = describePopulation(population)
n = numel(population);
if n == 0
    error('countModelOperations:EmptyPopulation', ...
        'population must contain at least one regressor.');
end
descriptors = repmat(factorizeGMPRegressor(population(1), 1), n, 1);
for index = 1:n
    descriptors(index) = factorizeGMPRegressor(population(index), index);
end
end

function dependencies = buildStandardDependencies(descriptors)
lags = zeros(1, 0);
powers = zeros(1, 0);
for index = 1:numel(descriptors)
    descriptor = descriptors(index);
    for term = 1:numel(descriptor.envelopeLags)
        [lags, powers] = ensurePower(lags, powers, ...
            descriptor.envelopeLags(term), ...
            descriptor.envelopePowers(term));
    end
end
[lags, order] = sort(lags);
dependencies.envelopeLags = lags;
dependencies.maximumPowers = powers(order);
dependencies.nEnvelopeLags = numel(lags);
end

function dependencies = buildStructureDependencies(descriptors, dependencies)
[dependencies.envelopeLags, dependencies.maximumPowers] = ...
    ensurePower(dependencies.envelopeLags, ...
    dependencies.maximumPowers, 0, 1);
for index = 1:numel(descriptors)
    descriptor = descriptors(index);
    if descriptor.canonicalGMP && descriptor.carrierLag == 0 && ...
            isempty(descriptor.envelopeLags)
        [dependencies.envelopeLags, dependencies.maximumPowers] = ...
            ensurePower(dependencies.envelopeLags, ...
            dependencies.maximumPowers, 0, 1);
    elseif descriptor.canonicalGMP && descriptor.carrierLag == 0 && ...
            all(descriptor.envelopeLags == 0)
        [dependencies.envelopeLags, dependencies.maximumPowers] = ...
            ensurePower(dependencies.envelopeLags, ...
            dependencies.maximumPowers, 0, ...
            sum(descriptor.envelopePowers) + 1);
    end
end
[dependencies.envelopeLags, order] = sort(dependencies.envelopeLags);
dependencies.maximumPowers = dependencies.maximumPowers(order);
canonical_mask = [descriptors.canonicalGMP];
carrier_lags = [descriptors(canonical_mask).carrierLag];
dependencies.nonzeroCarrierLags = unique( ...
    carrier_lags(carrier_lags ~= 0), 'sorted');
dependencies.nEnvelopeLags = numel(dependencies.envelopeLags);
end

function result = standardGenerationCost(descriptors, dependencies)
n_nonlinear = nnz(arrayfun(@(item) item.canonicalGMP && ...
    ~isempty(item.envelopeLags), descriptors));
result.M = 2*dependencies.nEnvelopeLags + ...
    sum(max(dependencies.maximumPowers - 1, 0)) + 2*n_nonlinear;
result.A = dependencies.nEnvelopeLags;
end

function result = structureGenerationCost(descriptors, dependencies)
formation = 0;
for index = 1:numel(descriptors)
    descriptor = descriptors(index);
    if ~descriptor.canonicalGMP || isempty(descriptor.envelopeLags)
        continue;
    end
    if descriptor.carrierLag == 0
        if any(descriptor.envelopeLags ~= 0)
            formation = formation + 1;
        end
    else
        formation = formation + 2;
    end
end
result.M = 2*dependencies.nEnvelopeLags + ...
    sum(max(dependencies.maximumPowers - 1, 0)) + formation;
result.A = dependencies.nEnvelopeLags;
end

function signatures = collectIQSignatures(descriptors)
signatures = strings(0, 1);
for index = 1:numel(descriptors)
    signatures(end+1, 1) = descriptors(index).iSignature; %#ok<AGROW>
    if ~descriptors(index).QColumnStructurallyZero
        signatures(end+1, 1) = descriptors(index).qSignature; %#ok<AGROW>
    end
end
end

function [lags, powers] = ensurePower(lags, powers, lag, power)
position = find(lags == lag, 1);
if isempty(position)
    lags(end+1) = lag;
    powers(end+1) = power;
else
    powers(position) = max(powers(position), power);
end
end

function support = validateSupport(support, population_size)
support = double(support(:));
if isempty(support) || any(~isfinite(support)) || ...
        any(support ~= floor(support)) || any(support < 1) || ...
        any(support > population_size) || ...
        numel(unique(support)) ~= numel(support)
    error('countModelOperations:InvalidSupport', ...
        'support must contain unique valid population indices.');
end
end
