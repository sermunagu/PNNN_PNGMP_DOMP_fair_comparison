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
assert(fixed.table.Properties.UserData.coefficientRangeDefinition == ...
    cfg.sweep.coefficientRangeDefinition);
assert(fixed.table.Properties.UserData.linearIdentificationScope == ...
    cfg.sweep.linearIdentificationScope);
assert(fixed.table.Properties.UserData.linearPrincipalLambda == 0);
assert(fixed.table.Properties.UserData.linearLambdaSelection == "none");
assert(fixed.table.Properties.UserData.fixedRidgeSupportPolicy == ...
    cfg.sweep.fixedRidgeSupportPolicy);
assert(isequal(fixed.paths, linear.paths));
assert(isequal(fixed.pnPathMap, linear.pnPathMap));
assert(height(fixed.table) == 18);
assert(width(fixed.table) == 8);
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
assert(all(isfinite(withPredictions.predictions.pnFull), 'all'));
assert(size(withPredictions.predictions.complexFull, 2) == ...
    numel(targets)*numel(cfg.fixedRidgeLambdas));
fixedSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'run_fixed_ridge_sweep.m'));
assert(~contains(fixedSource, 'selectDOMPSupport'));

%% Fixed Ridge ranges match explicit unit-peak, unit-column fits
[explicitComplex, explicitPN] = explicitFixedRanges( ...
    x, y, split, cfg, linear);
complexRows = string(withPredictions.table.Model) == "Complex GMP-DOMP";
pnRows = string(withPredictions.table.Model) == "PN-IQ PN-DOMP";
assert(all(abs(withPredictions.table.MaxAbsRealParameter(complexRows) - ...
    explicitComplex) < 1e-9));
assert(all(abs(withPredictions.table.MaxAbsRealParameter(pnRows) - ...
    explicitPN) < 1e-9));

%% Equivalent coefficient range and scientific metrics are scale invariant
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
assert(all(abs(scaled.predictions.pnFull - ...
    outputScale*withPredictions.predictions.pnFull) < 1e-8, 'all'));

fprintf('FIXED-LAMBDA LINEAR SWEEP TEST: PASS\n');

function [complexRanges, pnRanges] = explicitFixedRanges( ...
    x, y, split, cfg, linear)
rows = split.identificationIndices(:);
xNormalized = x / max(abs(x(rows)));
yNormalized = y / max(abs(y(rows)));
manager = GMP_createRegressorManager(xNormalized, yNormalized, cfg.gmp);
targets = cfg.sweep.parameterGrid(:);
lambdas = cfg.fixedRidgeLambdas(:);
variantCount = numel(targets)*numel(lambdas);
complexRanges = zeros(variantCount, 1);
pnRanges = zeros(variantCount, 1);

metadata = linear.pnPathMap(1:max(targets)/2, :);
complexSupport = unique(metadata.SourceRegressorIndex, 'stable');
complexRegressors = buildGMPRegressorRows( ...
    xNormalized, rows, manager, complexSupport);
rotation = complex(ones(numel(rows), 1));
nonzero = abs(xNormalized(rows)) ~= 0;
rotation(nonzero) = conj(xNormalized(rows(nonzero))) ./ ...
    abs(xNormalized(rows(nonzero)));
phaseNormalized = rotation .* complexRegressors;
[~, sourceColumns] = ismember(metadata.SourceRegressorIndex, complexSupport);
pnFeatures = zeros(numel(rows), height(metadata));
for featureIndex = 1:height(metadata)
    values = phaseNormalized(:, sourceColumns(featureIndex));
    if metadata.IsQ(featureIndex)
        pnFeatures(:, featureIndex) = imag(values);
    else
        pnFeatures(:, featureIndex) = real(values);
    end
end
rotatedTarget = rotation .* yNormalized(rows);

for targetIndex = 1:numel(targets)
    count = targets(targetIndex)/2;
    columns = (targetIndex - 1)*numel(lambdas) + (1:numel(lambdas));
    support = linear.paths.complex(1:count);
    regressors = buildGMPRegressorRows( ...
        xNormalized, rows, manager, support);
    regressors = regressors ./ vecnorm(regressors, 2, 1);
    features = pnFeatures(:, 1:count);
    features = features ./ vecnorm(features, 2, 1);
    for lambdaIndex = 1:numel(lambdas)
        lambda = lambdas(lambdaIndex);
        coefficients = ridgeFit(regressors, yNormalized(rows), lambda);
        complexRanges(columns(lambdaIndex)) = max([ ...
            abs(real(coefficients)); abs(imag(coefficients))]);
        coefficientsI = ridgeFit(features, real(rotatedTarget), lambda);
        coefficientsQ = ridgeFit(features, imag(rotatedTarget), lambda);
        pnRanges(columns(lambdaIndex)) = max(abs([ ...
            coefficientsI; coefficientsQ]));
    end
end
end

function coefficients = ridgeFit(regressors, target, lambda)
coefficients = (regressors'*regressors + ...
    lambda*eye(size(regressors, 2))) \ (regressors'*target);
end
