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
structuralZero = false(2*numel(population), 1);
for index = 1:numel(population)
    descriptor = factorizeGMPRegressor( ...
        manager.regPopulation(population(index)), population(index));
    structuralZero(numel(population) + index) = ...
        descriptor.QColumnStructurallyZero;
end
SourceRegressorIndex = [population; population];
IsQ = [false(numel(population), 1); true(numel(population), 1)];
pnFeatureMap = table(SourceRegressorIndex, IsQ);
pnFeatureMap = pnFeatureMap(~structuralZero, :);

complexPath = population(1:max(counts));
pnCandidates = find(~ismember( ...
    pnFeatureMap.SourceRegressorIndex, complexPath));
pnPath = pnCandidates(1:max(counts));
complexSupports = cell(numel(targets), 1);
pnFeatureSupports = cell(numel(targets), 1);
pnComplexSupports = cell(numel(targets), 1);
for index = 1:numel(targets)
    complexSupports{index} = complexPath(1:counts(index));
    pnFeatureSupports{index} = pnPath(1:counts(index));
    pnComplexSupports{index} = unique(pnFeatureMap.SourceRegressorIndex( ...
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
linear.paths = struct('complex', complexPath, 'pn', pnPath);
linear.pnPathMap = pnFeatureMap(pnPath, :);

fixed = run_fixed_ridge_sweep(x, y, split, cfg, linear);
assert(height(fixed.table) == 18);
assert(isequal(sort(unique(string(fixed.table.Model))), ...
    sort(["Complex GMP-DOMP"; "PN-IQ PN-DOMP"])));
assert(isequal(sort(unique(fixed.table.FixedLambda)), ...
    sort(cfg.fixedRidgeLambdas(:))));
for model = unique(string(fixed.table.Model)).'
    for lambda = cfg.fixedRidgeLambdas
        rows = string(fixed.table.Model) == model & ...
            fixed.table.FixedLambda == lambda;
        assert(nnz(rows) == numel(targets));
    end
end
assert(all(isfinite(fixed.table.IdentificationNMSEdB)));
assert(all(isfinite(fixed.table.FullSignalNMSEdB)));

withPredictions = run_fixed_ridge_sweep(x, y, split, cfg, linear, true);
assert(all(isfinite(withPredictions.predictions.complexFull), 'all'));
assert(all(isfinite(withPredictions.predictions.pnFull), 'all'));
assert(size(withPredictions.predictions.complexFull, 2) == ...
    numel(targets)*numel(cfg.fixedRidgeLambdas));

fprintf('FIXED-LAMBDA LINEAR SWEEP TEST: PASS\n');
