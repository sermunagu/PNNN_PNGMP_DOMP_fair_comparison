% Test fixed-ridge refits on stored, family-specific synthetic sweep supports.
% The fixture supplies paths directly and therefore invokes no DOMP or PNNN.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'metrics'));
addpath(fullfile(projectRoot, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

rng(731, 'twister');
n = 128;
x = 0.4*(randn(n, 1) + 1j*randn(n, 1));
y = 0.8*x + 0.12*x.*abs(x).^2 + ...
    0.005*(randn(n, 1) + 1j*randn(n, 1));
split.identificationIndices = (1:64).';
split.fullSignalIndices = (1:n).';
targets = [4; 6; 8];
counts = targets/2;

cfg = getFairDOMPComparisonConfig(projectRoot);
cfg.sweep.parameterGrid = targets.';
cfg.sweep.candidateBlockSize = 32;
cfg.gmp.blockSize = 32;
manager = GMP_createRegressorManager(x, y, cfg.gmp);
population = (1:numel(manager.regPopulation)).';
[probe, details] = buildPhaseNormalizedIQRegressors( ...
    x, 1, manager, population);
[~, reduction] = removeStructurallyZeroQFeatures( ...
    probe, details.featureMetadata, 1e-12);
featureMetadata = reduction.featureMetadata(reduction.keptIndices, :);

complexPath = population(1:max(counts));
pnCandidates = find(~ismember( ...
    featureMetadata.SourceRegressorIndex, complexPath));
pnPath = pnCandidates(1:max(counts));
complexSupports = cell(numel(targets), 1);
pnFeatureSupports = cell(numel(targets), 1);
pnComplexSupports = cell(numel(targets), 1);
for index = 1:numel(targets)
    complexSupports{index} = complexPath(1:counts(index));
    pnFeatureSupports{index} = pnPath(1:counts(index));
    pnComplexSupports{index} = unique(featureMetadata.SourceRegressorIndex( ...
        pnFeatureSupports{index}), 'stable');
end
assert(~isequal(complexSupports{end}, pnComplexSupports{end}));

TargetRealParameters = targets;
ActualRealParameters = targets;
FLOPsPerSample = 100 + targets;
linear.complexTable = table(TargetRealParameters, ...
    ActualRealParameters, FLOPsPerSample);
linear.pnTable = table(TargetRealParameters, ...
    ActualRealParameters, FLOPsPerSample);
linear.supports = struct('complex', {complexSupports}, ...
    'pnFeatures', {pnFeatureSupports}, 'pnComplex', {pnComplexSupports});
linear.paths = struct('complexTrain', flipud(complexPath), ...
    'complexIdentification', complexPath, 'pnTrain', flipud(pnPath), ...
    'pnIdentification', pnPath);
linear.featureMetadata = featureMetadata;

fixed = run_fixed_ridge_sweep(x, y, split, cfg, linear);
assert(height(fixed.table) == 18);
assert(isequal(sort(unique(string(fixed.table.Model))), ...
    sort(["Complex GMP-DOMP"; "PN-IQ PN-DOMP"])));
assert(isequal(sort(unique(fixed.table.FixedLambda)), ...
    [1e-5; 1e-4; 1e-3]));
for model = unique(string(fixed.table.Model)).'
    for lambda = fixed.fixedLambdas.'
        rows = string(fixed.table.Model) == model & ...
            fixed.table.FixedLambda == lambda;
        assert(nnz(rows) == numel(targets));
    end
end
assert(isequaln(fixed.supports.complex, complexSupports));
assert(isequaln(fixed.supports.pnFeatures, pnFeatureSupports));
assert(isequaln(fixed.supports.pnComplex, pnComplexSupports));
assert(isequaln(fixed.paths.complexIdentification, complexPath));
assert(isequaln(fixed.paths.pnIdentification, pnPath));
assert(fixed.metadata.dompInvocationCount == 0);
assert(fixed.metadata.pnnnTrainingCount == 0);
assert(all(structfun(@(value) value == 1, ...
    fixed.metadata.matrixPassCount)));
assert(all(isfinite(fixed.table.IdentificationNMSEdB)));
assert(all(isfinite(fixed.table.FullSignalNMSEdB)));

source = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'run_fixed_ridge_sweep.m'));
assert(~contains(source, 'selectDOMPSupport('));
assert(~contains(source, 'runPNNN'));
for removedWrapper = {'fitComplexVariants','fitPNVariants', ...
        'buildSelectedPNFeatures','evaluateComplexVariants', ...
        'evaluatePNVariants','buildResultTable'}
    assert(~contains(source, removedWrapper{1}));
end
assert(contains(source, 'conj(blockRotation) .* predictionNormalized'));

fprintf('FIXED-LAMBDA LINEAR SWEEP TEST: PASS\n');
