function fixed = run_fixed_ridge_sweep( ...
    x, y, split, cfg, linear, storePredictions)
% run_fixed_ridge_sweep - Refit stored linear supports with fixed ridge.
% The Complex GMP and PN-IQ families keep their own saved DOMP trajectories.
% Identification fits and one blockwise full-signal pass produce auxiliary rows.

if nargin < 6
    storePredictions = false;
end
if ~isscalar(storePredictions) || ~islogical(storePredictions)
    error('run_fixed_ridge_sweep:InvalidStorageMode', ...
        'storePredictions must be a logical scalar.');
end
fixedLambdas = [1e-3; 1e-4; 1e-5];
[targets, featureCounts, rows] = validateInputs(x, y, split, cfg, linear);
[referenceColumn, referenceTargetIndex] = referenceVariant( ...
    targets, fixedLambdas);
if storePredictions
    storedColumns = 1:(numel(targets)*numel(fixedLambdas));
else
    storedColumns = referenceColumn;
end
referenceStoredColumn = find(storedColumns == referenceColumn);
maximumFeatures = max(featureCounts);
complexPath = double(linear.paths.complexIdentification(:));
pnPath = double(linear.paths.pnIdentification(:));
complexSupport = complexPath(1:maximumFeatures);
pnFeaturePath = pnPath(1:maximumFeatures);
validateStoredSupports(linear, targets, featureCounts, ...
    complexPath, pnPath);

manager = GMP_createRegressorManager(x, y, cfg.gmp);
identificationTarget = y(rows.identification);
fullSignalTarget = y(rows.fullSignal);

fprintf('[Fixed ridge] Building one Complex GMP identification matrix...\n');
identificationU = buildGMPRegressorRows( ...
    x, rows.identification, manager, complexSupport);
[complexCoefficients, complexIdentificationNMSE, ...
    storedComplexIdentificationPredictions] = fitComplexVariants( ...
    identificationU, identificationTarget, featureCounts, fixedLambdas, ...
    storedColumns);
[complexFullNMSE, complexFullBuildCount, ...
    storedComplexFullPredictions] = ...
    evaluateComplexVariants( ...
    x, fullSignalTarget, rows.fullSignal, manager, complexSupport, ...
    complexCoefficients, cfg.gmp.blockSize, storedColumns);
clear identificationU complexCoefficients

fprintf('[Fixed ridge] Building one PN-IQ identification matrix...\n');
[identificationFeatures, identificationRotation, pnComplexSupport] = ...
    buildSelectedPNFeatures(x, rows.identification, manager, ...
    linear.featureMetadata, pnFeaturePath, cfg.sweep.candidateBlockSize);
[pnCoefficientsI, pnCoefficientsQ, pnIdentificationNMSE, ...
    storedPNIdentificationPredictions] = ...
    fitPNVariants(identificationFeatures, ...
    identificationRotation.*identificationTarget, identificationTarget, ...
    identificationRotation, featureCounts, fixedLambdas, storedColumns);
[pnFullNMSE, pnFullBuildCount, storedPNFullPredictions] = ...
    evaluatePNVariants( ...
    x, fullSignalTarget, rows.fullSignal, manager, ...
    linear.featureMetadata, pnFeaturePath, pnComplexSupport, ...
    pnCoefficientsI, pnCoefficientsQ, cfg.gmp.blockSize, storedColumns);

fixed.table = buildResultTable(linear, targets, fixedLambdas, ...
    complexIdentificationNMSE, complexFullNMSE, ...
    pnIdentificationNMSE, pnFullNMSE);
fixed.fixedLambdas = fixedLambdas;
fixed.supports = linear.supports;
fixed.paths = linear.paths;
referenceRow = fixed.table.Model == "Complex GMP-DOMP" & ...
    fixed.table.TargetRealParameters == targets(referenceTargetIndex) & ...
    fixed.table.FixedLambda == 1e-5;
fixed.reference = struct('model', "Complex GMP-DOMP", ...
    'targetRealParameters', targets(referenceTargetIndex), ...
    'fixedLambda', 1e-5, ...
    'support', linear.supports.complex{referenceTargetIndex}, ...
    'actualRealParameters', fixed.table.ActualRealParameters(referenceRow), ...
    'flopsPerSample', fixed.table.FLOPsPerSample(referenceRow), ...
    'identificationNMSEdB', ...
        fixed.table.IdentificationNMSEdB(referenceRow), ...
    'fullSignalNMSEdB', fixed.table.FullSignalNMSEdB(referenceRow), ...
    'identificationPrediction', ...
        storedComplexIdentificationPredictions(:, referenceStoredColumn), ...
    'fullSignalPrediction', ...
        storedComplexFullPredictions(:, referenceStoredColumn));
if storePredictions
    fixed.predictions = struct( ...
        'complexIdentification', storedComplexIdentificationPredictions, ...
        'complexFull', storedComplexFullPredictions, ...
        'pnIdentification', storedPNIdentificationPredictions, ...
        'pnFull', storedPNFullPredictions);
end
fixed.metadata = struct( ...
    'supportSource', "linear_sweep.mat", ...
    'supportContractsValidated', true, ...
    'dompInvocationCount', 0, ...
    'pnnnTrainingCount', 0, ...
    'matrixPassCount', struct('complexIdentification', 1, ...
        'complexFullSignal', 1, 'pnIdentification', 1, ...
        'pnFullSignal', 1), ...
    'fullSignalRegressorBuildCount', struct( ...
        'complex', complexFullBuildCount, 'pn', pnFullBuildCount), ...
    'fullSignalUsedForSelection', false, ...
    'fullSignalUsedForFitting', false);
end

function [column, targetIndex] = referenceVariant(targets, lambdas)
targetIndex = find(targets == 20, 1);
if isempty(targetIndex)
    targetIndex = 1;
end
lambdaIndex = find(lambdas == 1e-5, 1);
column = variantColumns(targetIndex, numel(lambdas));
column = column(lambdaIndex);
end

function [targets, featureCounts, rows] = validateInputs( ...
    x, y, split, cfg, linear)
x = x(:);
y = y(:);
requiredSplit = {'identificationIndices','fullSignalIndices'};
requiredLinear = {'complexTable','pnTable','supports','paths', ...
    'featureMetadata'};
if numel(x) ~= numel(y) || isempty(x) || any(~isfinite(x)) || ...
        any(~isfinite(y)) || ~all(isfield(split, requiredSplit)) || ...
        ~isstruct(linear) || ~all(isfield(linear, requiredLinear))
    error('run_fixed_ridge_sweep:InvalidInput', ...
        'Finite signals, evaluation rows, and the saved linear sweep are required.');
end
targets = double(cfg.sweep.parameterGrid(:));
if isempty(targets) || any(mod(targets, 2)) || ...
        ~isequal(targets, unique(targets, 'sorted'))
    error('run_fixed_ridge_sweep:InvalidGrid', ...
        'Targets must be sorted unique even parameter counts.');
end
featureCounts = targets/2;
rows = struct('identification', ...
    validateRows(split.identificationIndices, numel(x)), ...
    'fullSignal', validateRows(split.fullSignalIndices, numel(x)));
requiredMetadata = {'SourceRegressorIndex','Component'};
if ~istable(linear.featureMetadata) || ...
        ~all(ismember(requiredMetadata, ...
        linear.featureMetadata.Properties.VariableNames))
    error('run_fixed_ridge_sweep:InvalidMetadata', ...
        'The saved PN-IQ feature metadata is incomplete.');
end
end

function rows = validateRows(rows, upperBound)
rows = double(rows(:));
if isempty(rows) || any(~isfinite(rows)) || any(rows ~= floor(rows)) || ...
        any(rows < 1) || any(rows > upperBound)
    error('run_fixed_ridge_sweep:InvalidRows', ...
        'Split rows must be valid signal indices.');
end
end

function validateStoredSupports(linear, targets, counts, complexPath, pnPath)
if numel(complexPath) < max(counts) || numel(pnPath) < max(counts) || ...
        numel(linear.supports.complex) ~= numel(targets) || ...
        numel(linear.supports.pnFeatures) ~= numel(targets) || ...
        numel(linear.supports.pnComplex) ~= numel(targets)
    error('run_fixed_ridge_sweep:IncompleteSupports', ...
        'The saved paths and target-specific supports are incomplete.');
end
for index = 1:numel(targets)
    count = counts(index);
    expectedComplex = complexPath(1:count);
    expectedPNFeatures = pnPath(1:count);
    metadata = linear.featureMetadata(expectedPNFeatures, :);
    expectedPNComplex = unique(metadata.SourceRegressorIndex, 'stable');
    if ~isequal(double(linear.supports.complex{index}(:)), ...
            expectedComplex(:)) || ...
            ~isequal(double(linear.supports.pnFeatures{index}(:)), ...
            expectedPNFeatures(:)) || ...
            ~isequal(double(linear.supports.pnComplex{index}(:)), ...
            double(expectedPNComplex(:)))
        error('run_fixed_ridge_sweep:SupportMismatch', ...
            'A target support is not the prefix of its own saved trajectory.');
    end
end
end

function [coefficients, identificationNMSE, storedPredictions] = ...
    fitComplexVariants(U, target, counts, lambdas, storedColumns)
nVariants = numel(counts)*numel(lambdas);
coefficients = complex(zeros(max(counts), nVariants));
columnNorms = sqrt(sum(abs(U).^2, 1)).';
for targetIndex = 1:numel(counts)
    count = counts(targetIndex);
    columns = variantColumns(targetIndex, numel(lambdas));
    fit = fitComplexGMPGrid(U, target, 1:count, lambdas, columnNorms);
    coefficients(1:count, columns) = fit.coefficients;
end
predictions = U*coefficients;
identificationNMSE = nmseColumns(target, predictions);
storedPredictions = predictions(:, storedColumns);
end

function [coefficientsI, coefficientsQ, identificationNMSE, ...
    storedPredictions] = fitPNVariants( ...
    features, normalizedTarget, target, rotation, counts, lambdas, storedColumns)
featureNorms = sqrt(sum(features.^2, 1)).';
if any(featureNorms <= 0)
    error('run_fixed_ridge_sweep:ZeroPNFeature', ...
        'The stored PN-DOMP path contains a zero feature.');
end
normalizedFeatures = features ./ featureNorms.';
gram = normalizedFeatures.'*normalizedFeatures;
rhsI = normalizedFeatures.'*real(normalizedTarget);
rhsQ = normalizedFeatures.'*imag(normalizedTarget);
nVariants = numel(counts)*numel(lambdas);
coefficientsI = zeros(max(counts), nVariants);
coefficientsQ = coefficientsI;
for targetIndex = 1:numel(counts)
    count = counts(targetIndex);
    prefixGram = gram(1:count, 1:count);
    columns = variantColumns(targetIndex, numel(lambdas));
    for lambdaIndex = 1:numel(lambdas)
        regularizedGram = prefixGram + lambdas(lambdaIndex)*eye(count);
        coefficientsI(1:count, columns(lambdaIndex)) = ...
            (regularizedGram\rhsI(1:count))./featureNorms(1:count);
        coefficientsQ(1:count, columns(lambdaIndex)) = ...
            (regularizedGram\rhsQ(1:count))./featureNorms(1:count);
    end
end
prediction = conj(rotation).*(features*coefficientsI + ...
    1j*(features*coefficientsQ));
identificationNMSE = nmseColumns(target, prediction);
storedPredictions = prediction(:, storedColumns);
end

function columns = variantColumns(targetIndex, lambdaCount)
columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
end

function [features, rotation, supportComplex] = buildSelectedPNFeatures( ...
    x, rows, manager, metadata, path, blockSize)
selected = metadata(path, :);
supportComplex = unique(selected.SourceRegressorIndex, 'stable');
[found, columns] = ismember(selected.SourceRegressorIndex, supportComplex);
if any(~found)
    error('run_fixed_ridge_sweep:PNFeatureMapping', ...
        'The saved PN feature path cannot be mapped to its complex sources.');
end
isQ = string(selected.Component) == "Q";
columns(isQ) = columns(isQ) + numel(supportComplex);
features = zeros(numel(rows), numel(path));
for first = 1:blockSize:numel(rows)
    local = first:min(first + blockSize - 1, numel(rows));
    raw = buildPhaseNormalizedIQRegressors( ...
        x, rows(local), manager, supportComplex);
    features(local, :) = raw(:, columns);
end
rotation = computePhaseNormGMPRotation(x, rows);
end

function [nmse, buildCount, storedPredictions] = ...
    evaluateComplexVariants(x, target, rows, manager, support, ...
    coefficients, blockSize, storedColumns)
errorEnergy = zeros(1, size(coefficients, 2));
storedPredictions = complex(zeros(numel(rows), numel(storedColumns)));
buildCount = 0;
for first = 1:blockSize:numel(rows)
    local = first:min(first + blockSize - 1, numel(rows));
    U = buildGMPRegressorRows(x, rows(local), manager, support);
    prediction = U*coefficients;
    storedPredictions(local, :) = prediction(:, storedColumns);
    errorEnergy = errorEnergy + ...
        sum(abs(prediction - target(local)).^2, 1);
    buildCount = buildCount + 1;
end
nmse = energyToNMSE(errorEnergy, target);
end

function [nmse, buildCount, storedPredictions] = evaluatePNVariants( ...
    x, target, rows, manager, metadata, path, supportComplex, ...
    coefficientsI, coefficientsQ, blockSize, storedColumns)
selected = metadata(path, :);
[found, columns] = ismember(selected.SourceRegressorIndex, supportComplex);
if any(~found)
    error('run_fixed_ridge_sweep:PNFeatureMapping', ...
        'The full-signal PN feature path cannot be mapped.');
end
isQ = string(selected.Component) == "Q";
columns(isQ) = columns(isQ) + numel(supportComplex);
errorEnergy = zeros(1, size(coefficientsI, 2));
storedPredictions = complex(zeros(numel(rows), numel(storedColumns)));
buildCount = 0;
for first = 1:blockSize:numel(rows)
    local = first:min(first + blockSize - 1, numel(rows));
    [raw, details] = buildPhaseNormalizedIQRegressors( ...
        x, rows(local), manager, supportComplex);
    features = raw(:, columns);
    normalized = features*coefficientsI + 1j*(features*coefficientsQ);
    prediction = conj(details.phaseRotation).*normalized;
    storedPredictions(local, :) = prediction(:, storedColumns);
    errorEnergy = errorEnergy + ...
        sum(abs(prediction - target(local)).^2, 1);
    buildCount = buildCount + 1;
end
nmse = energyToNMSE(errorEnergy, target);
end

function values = nmseColumns(target, predictions)
values = energyToNMSE(sum(abs(predictions - target(:)).^2, 1), target);
end

function values = energyToNMSE(errorEnergy, target)
targetEnergy = sum(abs(target).^2);
if targetEnergy <= 0
    error('run_fixed_ridge_sweep:ZeroTargetEnergy', ...
        'NMSE requires nonzero target energy.');
end
values = 10*log10(errorEnergy(:)/targetEnergy);
end

function result = buildResultTable(linear, targets, lambdas, ...
    complexIdentificationNMSE, complexFullNMSE, ...
    pnIdentificationNMSE, pnFullNMSE)
lambdaCount = numel(lambdas);
variantTargets = repelem(targets(:), lambdaCount, 1);
variantLambdas = repmat(lambdas(:), numel(targets), 1);
complexActual = repeatPrincipalColumn( ...
    linear.complexTable, targets, 'ActualRealParameters', lambdaCount);
complexFLOPs = repeatPrincipalColumn( ...
    linear.complexTable, targets, 'FLOPsPerSample', lambdaCount);
pnActual = repeatPrincipalColumn( ...
    linear.pnTable, targets, 'ActualRealParameters', lambdaCount);
pnFLOPs = repeatPrincipalColumn( ...
    linear.pnTable, targets, 'FLOPsPerSample', lambdaCount);
Model = [repmat("Complex GMP-DOMP", numel(variantTargets), 1); ...
    repmat("PN-IQ PN-DOMP", numel(variantTargets), 1)];
TargetRealParameters = [variantTargets; variantTargets];
ActualRealParameters = [complexActual; pnActual];
FixedLambda = [variantLambdas; variantLambdas];
IdentificationNMSEdB = [complexIdentificationNMSE; pnIdentificationNMSE];
FullSignalNMSEdB = [complexFullNMSE; pnFullNMSE];
FLOPsPerSample = [complexFLOPs; pnFLOPs];
result = table(Model, TargetRealParameters, ActualRealParameters, ...
    FixedLambda, IdentificationNMSEdB, FullSignalNMSEdB, FLOPsPerSample);
end

function values = repeatPrincipalColumn(principal, targets, name, repeatCount)
values = zeros(numel(targets), 1);
for index = 1:numel(targets)
    row = principal.TargetRealParameters == targets(index);
    if nnz(row) ~= 1
        error('run_fixed_ridge_sweep:PrincipalRowMismatch', ...
            'Each target must have exactly one principal-family row.');
    end
    values(index) = principal.(name)(row);
end
values = repelem(values, repeatCount, 1);
end
