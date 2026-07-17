function fixed = run_fixed_ridge_sweep( ...
    x, y, split, cfg, linear, storePredictions)
% run_fixed_ridge_sweep - Refit saved linear supports with fixed Ridge.
% Complex GMP and PN-IQ retain their own DOMP paths; Ridge changes only
% coefficient estimation and never invokes DOMP or trains a PNNN.

if nargin < 6
    storePredictions = false;
end
x = x(:);
y = y(:);
targets = double(cfg.sweep.parameterGrid(:));
featureCounts = targets/2;
maximumFeatures = max(featureCounts);
identificationRows = double(split.identificationIndices(:));
fullSignalRows = double(split.fullSignalIndices(:));

fixedLambdas = cfg.fixedRidgeLambdas(:);
lambdaCount = numel(fixedLambdas);
variantCount = numel(targets)*lambdaCount;
complexPath = double(linear.paths.complexIdentification(:));
pnPath = double(linear.paths.pnIdentification(:));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    expectedComplex = complexPath(1:count);
    expectedPNFeatures = pnPath(1:count);
    metadata = linear.featureMetadata(expectedPNFeatures, :);
    expectedPNComplex = unique(metadata.SourceRegressorIndex, 'stable');
    if ~isequal(double(linear.supports.complex{targetIndex}(:)), ...
            expectedComplex) || ...
            ~isequal(double(linear.supports.pnFeatures{targetIndex}(:)), ...
            expectedPNFeatures) || ...
            ~isequal(double(linear.supports.pnComplex{targetIndex}(:)), ...
            double(expectedPNComplex(:)))
        error('run_fixed_ridge_sweep:SupportMismatch', ...
            'Every target must be the exact prefix of its own saved path.');
    end
end

referenceTargetIndex = 1;
referenceColumn = (referenceTargetIndex - 1)*lambdaCount + lambdaCount;
if storePredictions
    storedColumns = 1:variantCount;
else
    storedColumns = referenceColumn;
end
referenceStoredColumn = find(storedColumns == referenceColumn, 1);

manager = GMP_createRegressorManager(x, y, cfg.gmp);
identificationTarget = y(identificationRows);
fullSignalTarget = y(fullSignalRows);
targetEnergyIdentification = sum(abs(identificationTarget).^2);
targetEnergyFull = sum(abs(fullSignalTarget).^2);

%% Complex GMP: normalize, fit all fixed lambdas, and predict
% The complex column norm is computed once on identification and undone in h.
fprintf('[Fixed ridge] Building one Complex GMP identification matrix...\n');
complexSupport = complexPath(1:maximumFeatures);
identificationU = buildGMPRegressorRows( ...
    x, identificationRows, manager, complexSupport);
complexColumnNorms = sqrt(sum(abs(identificationU).^2, 1)).';
complexCoefficients = complex(zeros(maximumFeatures, variantCount));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    columnNorms = complexColumnNorms(1:count);
    normalizedU = identificationU(:, 1:count) ./ columnNorms.';
    gram = normalizedU' * normalizedU;
    rhs = normalizedU' * identificationTarget;
    for lambdaIndex = 1:lambdaCount
        normalizedCoefficients = ...
            (gram + fixedLambdas(lambdaIndex)*eye(count)) \ rhs;
        complexCoefficients(1:count, columns(lambdaIndex)) = ...
            normalizedCoefficients ./ columnNorms;
    end
end
complexIdentificationPredictions = identificationU * complexCoefficients;
complexIdentificationError = sum(abs( ...
    complexIdentificationPredictions - identificationTarget).^2, 1);
complexIdentificationNMSE = ...
    10*log10(complexIdentificationError(:)/targetEnergyIdentification);
storedComplexIdentificationPredictions = ...
    complexIdentificationPredictions(:, storedColumns);

complexFullError = zeros(1, variantCount);
storedComplexFullPredictions = ...
    complex(zeros(numel(fullSignalRows), numel(storedColumns)));
complexFullBuildCount = 0;
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    U = buildGMPRegressorRows( ...
        x, fullSignalRows(local), manager, complexSupport);
    prediction = U * complexCoefficients;
    storedComplexFullPredictions(local, :) = prediction(:, storedColumns);
    complexFullError = complexFullError + ...
        sum(abs(prediction - fullSignalTarget(local)).^2, 1);
    complexFullBuildCount = complexFullBuildCount + 1;
end
complexFullNMSE = 10*log10(complexFullError(:)/targetEnergyFull);
clear identificationU complexCoefficients complexIdentificationPredictions

%% PN-IQ: normalize real features, fit I/Q, and restore phase
% Both outputs use the same real features; only their coefficients differ.
fprintf('[Fixed ridge] Building one PN-IQ identification matrix...\n');
pnFeaturePath = pnPath(1:maximumFeatures);
selectedMetadata = linear.featureMetadata(pnFeaturePath, :);
pnComplexSupport = unique(selectedMetadata.SourceRegressorIndex, 'stable');
[~, selectedColumns] = ismember( ...
    selectedMetadata.SourceRegressorIndex, pnComplexSupport);
isQ = string(selectedMetadata.Component) == "Q";
selectedColumns(isQ) = selectedColumns(isQ) + numel(pnComplexSupport);
pnDescriptors = repmat(factorizeGMPRegressor( ...
    manager.regPopulation(pnComplexSupport(1)), pnComplexSupport(1)), ...
    numel(pnComplexSupport), 1);
for index = 1:numel(pnComplexSupport)
    pnDescriptors(index) = factorizeGMPRegressor( ...
        manager.regPopulation(pnComplexSupport(index)), ...
        pnComplexSupport(index));
end
identificationInput = x(identificationRows);
identificationRotation = complex(ones(numel(identificationRows), 1));
nonzero = abs(identificationInput) ~= 0;
identificationRotation(nonzero) = ...
    conj(identificationInput(nonzero)) ./ abs(identificationInput(nonzero));

identificationFeatures = zeros(numel(identificationRows), maximumFeatures);
for first = 1:cfg.sweep.candidateBlockSize:numel(identificationRows)
    local = first:min(first + cfg.sweep.candidateBlockSize - 1, ...
        numel(identificationRows));
    blockRows = identificationRows(local);
    blockRotation = identificationRotation(local);
    complexRegressors = buildGMPRegressorRows( ...
        x, blockRows, manager, pnComplexSupport);
    phaseNormalized = blockRotation .* complexRegressors;
    regressorsI = zeros(numel(local), numel(pnComplexSupport));
    regressorsQ = zeros(numel(local), numel(pnComplexSupport));
    for regressorIndex = 1:numel(pnComplexSupport)
        descriptor = pnDescriptors(regressorIndex);
        if descriptor.canonicalGMP
            carrierRows = mod(blockRows - descriptor.carrierLag - 1, ...
                numel(x)) + 1;
            normalizedCarrier = blockRotation .* x(carrierRows);
            envelope = ones(numel(local), 1);
            for termIndex = 1:numel(descriptor.envelopeLags)
                envelopeRows = mod(blockRows - ...
                    descriptor.envelopeLags(termIndex) - 1, numel(x)) + 1;
                envelope = envelope .* abs(x(envelopeRows)).^ ...
                    descriptor.envelopePowers(termIndex);
            end
            regressorsI(:, regressorIndex) = ...
                real(normalizedCarrier) .* envelope;
            regressorsQ(:, regressorIndex) = ...
                imag(normalizedCarrier) .* envelope;
            if descriptor.QColumnStructurallyZero
                regressorsQ(:, regressorIndex) = 0;
            end
        else
            regressorsI(:, regressorIndex) = ...
                real(phaseNormalized(:, regressorIndex));
            regressorsQ(:, regressorIndex) = ...
                imag(phaseNormalized(:, regressorIndex));
        end
    end
    raw = [regressorsI, regressorsQ];
    identificationFeatures(local, :) = raw(:, selectedColumns);
end
normalizedTarget = identificationRotation .* identificationTarget;

featureNorms = sqrt(sum(identificationFeatures.^2, 1)).';
normalizedFeatures = identificationFeatures ./ featureNorms.';
gram = normalizedFeatures.' * normalizedFeatures;
rhsI = normalizedFeatures.' * real(normalizedTarget);
rhsQ = normalizedFeatures.' * imag(normalizedTarget);
pnCoefficientsI = zeros(maximumFeatures, variantCount);
pnCoefficientsQ = zeros(maximumFeatures, variantCount);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    prefixGram = gram(1:count, 1:count);
    for lambdaIndex = 1:lambdaCount
        regularizedGram = ...
            prefixGram + fixedLambdas(lambdaIndex)*eye(count);
        pnCoefficientsI(1:count, columns(lambdaIndex)) = ...
            (regularizedGram \ rhsI(1:count)) ./ featureNorms(1:count);
        pnCoefficientsQ(1:count, columns(lambdaIndex)) = ...
            (regularizedGram \ rhsQ(1:count)) ./ featureNorms(1:count);
    end
end
pnIdentificationNormalized = ...
    identificationFeatures * pnCoefficientsI + ...
    1j * (identificationFeatures * pnCoefficientsQ);
pnIdentificationPredictions = ...
    conj(identificationRotation) .* pnIdentificationNormalized;
pnIdentificationError = sum(abs( ...
    pnIdentificationPredictions - identificationTarget).^2, 1);
pnIdentificationNMSE = ...
    10*log10(pnIdentificationError(:)/targetEnergyIdentification);
storedPNIdentificationPredictions = ...
    pnIdentificationPredictions(:, storedColumns);

pnFullError = zeros(1, variantCount);
storedPNFullPredictions = ...
    complex(zeros(numel(fullSignalRows), numel(storedColumns)));
pnFullBuildCount = 0;
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    blockRows = fullSignalRows(local);
    blockInput = x(blockRows);
    blockRotation = complex(ones(numel(local), 1));
    nonzero = abs(blockInput) ~= 0;
    blockRotation(nonzero) = ...
        conj(blockInput(nonzero)) ./ abs(blockInput(nonzero));
    complexRegressors = buildGMPRegressorRows( ...
        x, blockRows, manager, pnComplexSupport);
    phaseNormalized = blockRotation .* complexRegressors;
    regressorsI = zeros(numel(local), numel(pnComplexSupport));
    regressorsQ = zeros(numel(local), numel(pnComplexSupport));
    for regressorIndex = 1:numel(pnComplexSupport)
        descriptor = pnDescriptors(regressorIndex);
        if descriptor.canonicalGMP
            carrierRows = mod(blockRows - descriptor.carrierLag - 1, ...
                numel(x)) + 1;
            normalizedCarrier = blockRotation .* x(carrierRows);
            envelope = ones(numel(local), 1);
            for termIndex = 1:numel(descriptor.envelopeLags)
                envelopeRows = mod(blockRows - ...
                    descriptor.envelopeLags(termIndex) - 1, numel(x)) + 1;
                envelope = envelope .* abs(x(envelopeRows)).^ ...
                    descriptor.envelopePowers(termIndex);
            end
            regressorsI(:, regressorIndex) = ...
                real(normalizedCarrier) .* envelope;
            regressorsQ(:, regressorIndex) = ...
                imag(normalizedCarrier) .* envelope;
            if descriptor.QColumnStructurallyZero
                regressorsQ(:, regressorIndex) = 0;
            end
        else
            regressorsI(:, regressorIndex) = ...
                real(phaseNormalized(:, regressorIndex));
            regressorsQ(:, regressorIndex) = ...
                imag(phaseNormalized(:, regressorIndex));
        end
    end
    raw = [regressorsI, regressorsQ];
    features = raw(:, selectedColumns);
    predictionNormalized = complex( ...
        features * pnCoefficientsI, features * pnCoefficientsQ);
    prediction = conj(blockRotation) .* predictionNormalized;
    storedPNFullPredictions(local, :) = prediction(:, storedColumns);
    pnFullError = pnFullError + ...
        sum(abs(prediction - fullSignalTarget(local)).^2, 1);
    pnFullBuildCount = pnFullBuildCount + 1;
end
pnFullNMSE = 10*log10(pnFullError(:)/targetEnergyFull);

%% Package the two supplementary Ridge families
% Parameters and FLOPs are copied from the principal rows because support is fixed.
variantTargets = repelem(targets, lambdaCount, 1);
variantLambdas = repmat(fixedLambdas, numel(targets), 1);
complexActual = zeros(numel(targets), 1);
complexFLOPs = zeros(numel(targets), 1);
pnActual = zeros(numel(targets), 1);
pnFLOPs = zeros(numel(targets), 1);
for targetIndex = 1:numel(targets)
    complexRow = linear.complexTable.TargetRealParameters == targets(targetIndex);
    pnRow = linear.pnTable.TargetRealParameters == targets(targetIndex);
    complexActual(targetIndex) = ...
        linear.complexTable.ActualRealParameters(complexRow);
    complexFLOPs(targetIndex) = linear.complexTable.FLOPsPerSample(complexRow);
    pnActual(targetIndex) = linear.pnTable.ActualRealParameters(pnRow);
    pnFLOPs(targetIndex) = linear.pnTable.FLOPsPerSample(pnRow);
end
complexActual = repelem(complexActual, lambdaCount, 1);
complexFLOPs = repelem(complexFLOPs, lambdaCount, 1);
pnActual = repelem(pnActual, lambdaCount, 1);
pnFLOPs = repelem(pnFLOPs, lambdaCount, 1);

Model = [repmat("Complex GMP-DOMP", variantCount, 1); ...
    repmat("PN-IQ PN-DOMP", variantCount, 1)];
TargetRealParameters = [variantTargets; variantTargets];
ActualRealParameters = [complexActual; pnActual];
FixedLambda = [variantLambdas; variantLambdas];
IdentificationNMSEdB = [complexIdentificationNMSE; pnIdentificationNMSE];
FullSignalNMSEdB = [complexFullNMSE; pnFullNMSE];
FLOPsPerSample = [complexFLOPs; pnFLOPs];
fixed.table = table(Model, TargetRealParameters, ActualRealParameters, ...
    FixedLambda, IdentificationNMSEdB, FullSignalNMSEdB, FLOPsPerSample);
fixed.fixedLambdas = fixedLambdas;
fixed.supports = linear.supports;
fixed.paths = linear.paths;

referenceRow = fixed.table.Model == "Complex GMP-DOMP" & ...
    fixed.table.TargetRealParameters == targets(referenceTargetIndex) & ...
    fixed.table.FixedLambda == fixedLambdas(end);
fixed.reference = struct('model', "Complex GMP-DOMP", ...
    'targetRealParameters', targets(referenceTargetIndex), ...
    'fixedLambda', fixedLambdas(end), ...
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
fixed.metadata = struct('supportSource', "linear_sweep.mat", ...
    'supportContractsValidated', true, 'dompInvocationCount', 0, ...
    'pnnnTrainingCount', 0, ...
    'matrixPassCount', struct('complexIdentification', 1, ...
    'complexFullSignal', 1, 'pnIdentification', 1, 'pnFullSignal', 1), ...
    'fullSignalRegressorBuildCount', struct( ...
    'complex', complexFullBuildCount, 'pn', pnFullBuildCount), ...
    'fullSignalUsedForSelection', false, ...
    'fullSignalUsedForFitting', false);
end
