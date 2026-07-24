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
pniqFeatureMap = table(SourceRegressorIndex, IsQ);
pniqFeatureMap = pniqFeatureMap(~structuralZero, :);

complexPath = population(1:max(counts));
pniqCandidates = find(~ismember( ...
    pniqFeatureMap.SourceRegressorIndex, complexPath));
pniqPath = pniqCandidates(1:max(counts));
complexSupports = cell(numel(targets), 1);
pniqFeatureSupports = cell(numel(targets), 1);
pniqComplexSupports = cell(numel(targets), 1);
for index = 1:numel(targets)
    complexSupports{index} = complexPath(1:counts(index));
    pniqFeatureSupports{index} = pniqPath(1:counts(index));
    pniqComplexSupports{index} = unique( ...
        pniqFeatureMap.SourceRegressorIndex( ...
        pniqFeatureSupports{index}), 'stable');
end
assert(~isequal(complexSupports{end}, pniqComplexSupports{end}));

TargetRealParameters = targets;
ActualRealParameters = targets;
FLOPsPerSample = 100 + targets;
linear.complexTable = table(TargetRealParameters, ...
    ActualRealParameters, FLOPsPerSample);
linear.pniqTable = table(TargetRealParameters, ...
    ActualRealParameters, FLOPsPerSample);
linear.paths = struct('complex', complexPath, 'pniq', pniqPath);
linear.pniqPathMap = pniqFeatureMap(pniqPath, :);

fixed = run_fixed_ridge_sweep(x, y, split, cfg, linear);
assert(fixed.table.Properties.UserData.coefficientRangeDefinition == ...
    cfg.sweep.coefficientRangeDefinition);
assert(fixed.table.Properties.UserData.linearIdentificationScope == ...
    cfg.sweep.linearIdentificationScope);
assert(fixed.table.Properties.UserData.linearPrincipalLambda == 0);
assert(fixed.table.Properties.UserData.linearLambdaSelection == "none");
assert(fixed.table.Properties.UserData.fixedRidgeSupportPolicy == ...
    cfg.sweep.fixedRidgeSupportPolicy);
assert(isequal(fixed.paths, linear.paths));
assert(isequal(fixed.pniqPathMap, linear.pniqPathMap));
assert(height(fixed.table) == 18);
assert(width(fixed.table) == 8);
assert(isequal(sort(unique(string(fixed.table.Model))), ...
    sort([cfg.names.complexGMPDOMP; cfg.names.pniqGMP])));
assert(isequal(sort(unique(fixed.table.FixedLambda)), ...
    sort(cfg.fixedRidgeLambdas(:))));
for model = unique(string(fixed.table.Model)).'
    for lambda = cfg.fixedRidgeLambdas
        rows = string(fixed.table.Model) == model & ...
            fixed.table.FixedLambda == lambda;
        assert(nnz(rows) == numel(targets));
    end
    for target = targets.'
        rows = string(fixed.table.Model) == model & ...
            fixed.table.TargetRealParameters == target;
        assert(isscalar(unique(fixed.table.ActualRealParameters(rows))));
        assert(isscalar(unique(fixed.table.FLOPsPerSample(rows))));
    end
end
assert(all(isfinite(fixed.table.IdentificationNMSEdB)));
assert(all(isfinite(fixed.table.FullSignalNMSEdB)));
assert(all(isfinite(fixed.table.MaxAbsRealParameter)));
assert(all(fixed.table.MaxAbsRealParameter >= 0));

withPredictions = run_fixed_ridge_sweep(x, y, split, cfg, linear, true);
assert(all(isfinite(withPredictions.predictions.complexFull), 'all'));
assert(all(isfinite(withPredictions.predictions.pniqFull), 'all'));
assert(size(withPredictions.predictions.complexFull, 2) == ...
    numel(targets)*numel(cfg.fixedRidgeLambdas));
fixedSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'run_fixed_ridge_sweep.m'));
assert(~contains(fixedSource, 'selectDOMPSupport'));

%% Fixed Ridge ranges match explicit per-column peak fits
[explicitComplex, explicitPNIQ] = explicitFixedRanges( ...
    x, y, split, cfg, linear);
complexRows = string(withPredictions.table.Model) == cfg.names.complexGMPDOMP;
pniqRows = string(withPredictions.table.Model) == cfg.names.pniqGMP;
assert(all(abs(withPredictions.table.MaxAbsRealParameter(complexRows) - ...
    explicitComplex) < 1e-9));
assert(all(abs(withPredictions.table.MaxAbsRealParameter(pniqRows) - ...
    explicitPNIQ) < 1e-9));

%% Equivalent normalized coefficients are invariant to input/output scaling
inputScale = 1.7;
outputScale = 0.6;
scaled = run_fixed_ridge_sweep(inputScale*x, outputScale*y, split, ...
    cfg, linear, true);
assert(isequal(scaled.table(:, {'Model','TargetRealParameters', ...
    'ActualRealParameters','FixedLambda','FLOPsPerSample'}), ...
    withPredictions.table(:, {'Model','TargetRealParameters', ...
    'ActualRealParameters','FixedLambda','FLOPsPerSample'})));
nmseColumns = {'IdentificationNMSEdB','FullSignalNMSEdB'};
originalNMSE = withPredictions.table{:, nmseColumns};
scaledNMSE = scaled.table{:, nmseColumns};
assert(all(abs(scaledNMSE - originalNMSE) < 1e-8 | ...
    (scaledNMSE < -250 & originalNMSE < -250), 'all'));
assert(all(abs(scaled.table.MaxAbsRealParameter - ...
    withPredictions.table.MaxAbsRealParameter) < 1e-8));
assert(all(abs(scaled.predictions.complexFull - ...
    outputScale*withPredictions.predictions.complexFull) < 1e-8, 'all'));
assert(all(abs(scaled.predictions.pniqFull - ...
    outputScale*withPredictions.predictions.pniqFull) < 1e-8, 'all'));

fprintf('FIXED-LAMBDA LINEAR SWEEP TEST: PASS\n');

function [complexRanges, pniqRanges] = explicitFixedRanges( ...
    x, y, split, cfg, linear)
rows = split.identificationIndices(:);
outputPeak = max(abs(y(rows)));
unitPeakTarget = y / outputPeak;
manager = GMP_createRegressorManager(x, y, cfg.gmp);
targets = cfg.sweep.parameterGrid(:);
lambdas = cfg.fixedRidgeLambdas(:);
variantCount = numel(targets)*numel(lambdas);
complexRanges = zeros(variantCount, 1);
pniqRanges = zeros(variantCount, 1);

metadata = linear.pniqPathMap(1:max(targets)/2, :);
complexSupport = unique(metadata.SourceRegressorIndex, 'stable');
complexRegressors = buildGMPRegressorRows( ...
    x, rows, manager, complexSupport);
rotation = complex(ones(numel(rows), 1));
nonzero = abs(x(rows)) ~= 0;
rotation(nonzero) = conj(x(rows(nonzero))) ./ abs(x(rows(nonzero)));
phaseNormalized = rotation .* complexRegressors;
[~, sourceColumns] = ismember(metadata.SourceRegressorIndex, complexSupport);
pniqFeatures = zeros(numel(rows), height(metadata));
for featureIndex = 1:height(metadata)
    values = phaseNormalized(:, sourceColumns(featureIndex));
    if metadata.IsQ(featureIndex)
        pniqFeatures(:, featureIndex) = imag(values);
    else
        pniqFeatures(:, featureIndex) = real(values);
    end
end
phaseNormalizedTarget = rotation .* unitPeakTarget(rows);

for targetIndex = 1:numel(targets)
    count = targets(targetIndex)/2;
    columns = (targetIndex - 1)*numel(lambdas) + (1:numel(lambdas));
    support = linear.paths.complex(1:count);
    regressors = buildGMPRegressorRows( ...
        x, rows, manager, support);
    regressors = regressors ./ max(abs(regressors), [], 1);
    features = pniqFeatures(:, 1:count);
    features = features ./ max(abs(features), [], 1);
    for lambdaIndex = 1:numel(lambdas)
        lambda = lambdas(lambdaIndex);
        coefficients = ridgeFit(regressors, unitPeakTarget(rows), lambda);
        complexRanges(columns(lambdaIndex)) = max(abs(coefficients));
        coefficientsI = ridgeFit(features, real(phaseNormalizedTarget), lambda);
        coefficientsQ = ridgeFit(features, imag(phaseNormalizedTarget), lambda);
        pniqRanges(columns(lambdaIndex)) = max(abs([ ...
            coefficientsI; coefficientsQ]));
    end
end
end

function coefficients = ridgeFit(regressors, target, lambda)
coefficients = (regressors'*regressors + ...
    lambda*eye(size(regressors, 2))) \ (regressors'*target);
end
