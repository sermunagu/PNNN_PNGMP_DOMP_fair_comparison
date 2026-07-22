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
complexPath = double(linear.paths.complex(:));
if storePredictions
    storedColumns = 1:variantCount;
else
    storedColumns = zeros(1, 0);
end

manager = GMP_createRegressorManager(x, y, cfg.gmp);
identificationTarget = y(identificationRows);
fullSignalTarget = y(fullSignalRows);
targetEnergyIdentification = sum(abs(identificationTarget).^2);
targetEnergyFull = sum(abs(fullSignalTarget).^2);
outputPeak = max(abs(y(identificationRows)));

%% Complex GMP: normalize, fit all fixed lambdas, and predict
% The complex column norm is computed once on identification and undone in h.
% Unit-peak input gives the same unit-norm columns because its global scale
% cancels when each homogeneous GMP column is normalized.
fprintf('[Fixed ridge] Building one %s identification matrix...\n', ...
    cfg.names.complexGMPDOMP);
complexSupport = complexPath(1:maximumFeatures);
identificationU = buildGMPRegressorRows( ...
    x, identificationRows, manager, complexSupport);
complexColumnNorms = sqrt(sum(abs(identificationU).^2, 1)).';
complexCoefficients = complex(zeros(maximumFeatures, variantCount));
complexComparisonCoefficients = complexCoefficients;
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
        complexComparisonCoefficients(1:count, columns(lambdaIndex)) = ...
            normalizedCoefficients / outputPeak;
        complexCoefficients(1:count, columns(lambdaIndex)) = ...
            normalizedCoefficients ./ columnNorms;
    end
end
complexIdentificationPredictions = identificationU * complexCoefficients;
complexIdentificationError = sum(abs( ...
    complexIdentificationPredictions - identificationTarget).^2, 1);
complexIdentificationNMSE = ...
    10*log10(complexIdentificationError(:)/targetEnergyIdentification);
complexMaxAbs = zeros(variantCount, 1);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    for lambdaIndex = 1:lambdaCount
        values = complexComparisonCoefficients( ...
            1:count, columns(lambdaIndex));
        complexMaxAbs(columns(lambdaIndex)) = max([ ...
            abs(real(values)); abs(imag(values))]);
    end
end

complexFullError = zeros(1, variantCount);
storedComplexFullPredictions = complex(zeros(numel(fullSignalRows), numel(storedColumns)));
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    U = buildGMPRegressorRows(x, fullSignalRows(local), manager, complexSupport);
    prediction = U * complexCoefficients;
    storedComplexFullPredictions(local, :) = prediction(:, storedColumns);
    complexFullError = complexFullError + sum(abs(prediction - fullSignalTarget(local)).^2, 1);
end
complexFullNMSE = 10*log10(complexFullError(:)/targetEnergyFull);
clear identificationU complexCoefficients complexIdentificationPredictions

%% PN-IQ: normalize real features, fit I/Q, and restore phase
% Both outputs use the same real features; only their coefficients differ.
fprintf('[Fixed ridge] Building one %s identification matrix...\n', ...
    cfg.names.pniqGMP);
selectedMetadata = linear.pniqPathMap(1:maximumFeatures, :);
pniqComplexSupport = unique(selectedMetadata.SourceRegressorIndex, 'stable');
[~, selectedColumns] = ismember( ...
    selectedMetadata.SourceRegressorIndex, pniqComplexSupport);
isQ = selectedMetadata.IsQ;
selectedColumns(isQ) = selectedColumns(isQ) + numel(pniqComplexSupport);
pniqDescriptors = repmat(factorizeGMPRegressor( ...
    manager.regPopulation(pniqComplexSupport(1)), pniqComplexSupport(1)), ...
    numel(pniqComplexSupport), 1);
for index = 1:numel(pniqComplexSupport)
    pniqDescriptors(index) = factorizeGMPRegressor( ...
        manager.regPopulation(pniqComplexSupport(index)), ...
        pniqComplexSupport(index));
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
        x, blockRows, manager, pniqComplexSupport);
    phaseNormalized = blockRotation .* complexRegressors;
    regressorsI = zeros(numel(local), numel(pniqComplexSupport));
    regressorsQ = zeros(numel(local), numel(pniqComplexSupport));
    for regressorIndex = 1:numel(pniqComplexSupport)
        descriptor = pniqDescriptors(regressorIndex);
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
pniqCoefficientsI = zeros(maximumFeatures, variantCount);
pniqCoefficientsQ = zeros(maximumFeatures, variantCount);
pniqComparisonI = pniqCoefficientsI;
pniqComparisonQ = pniqCoefficientsQ;
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    prefixGram = gram(1:count, 1:count);
    for lambdaIndex = 1:lambdaCount
        regularizedGram = ...
            prefixGram + fixedLambdas(lambdaIndex)*eye(count);
        normalizedI = regularizedGram \ rhsI(1:count);
        normalizedQ = regularizedGram \ rhsQ(1:count);
        pniqCoefficientsI(1:count, columns(lambdaIndex)) = ...
            normalizedI ./ featureNorms(1:count);
        pniqCoefficientsQ(1:count, columns(lambdaIndex)) = ...
            normalizedQ ./ featureNorms(1:count);
        pniqComparisonI(1:count, columns(lambdaIndex)) = ...
            normalizedI / outputPeak;
        pniqComparisonQ(1:count, columns(lambdaIndex)) = ...
            normalizedQ / outputPeak;
    end
end
pniqIdentificationNormalized = ...
    identificationFeatures * pniqCoefficientsI + ...
    1j * (identificationFeatures * pniqCoefficientsQ);
pniqIdentificationPredictions = ...
    conj(identificationRotation) .* pniqIdentificationNormalized;
pniqIdentificationError = sum(abs( ...
    pniqIdentificationPredictions - identificationTarget).^2, 1);
pniqIdentificationNMSE = ...
    10*log10(pniqIdentificationError(:)/targetEnergyIdentification);
pniqMaxAbs = zeros(variantCount, 1);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    for lambdaIndex = 1:lambdaCount
        column = columns(lambdaIndex);
        pniqMaxAbs(column) = max(abs([ ...
            pniqComparisonI(1:count, column); ...
            pniqComparisonQ(1:count, column)]));
    end
end

pniqFullError = zeros(1, variantCount);
storedPNIQFullPredictions = ...
    complex(zeros(numel(fullSignalRows), numel(storedColumns)));
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    blockRows = fullSignalRows(local);
    blockInput = x(blockRows);
    blockRotation = complex(ones(numel(local), 1));
    nonzero = abs(blockInput) ~= 0;
    blockRotation(nonzero) = ...
        conj(blockInput(nonzero)) ./ abs(blockInput(nonzero));
    complexRegressors = buildGMPRegressorRows( ...
        x, blockRows, manager, pniqComplexSupport);
    phaseNormalized = blockRotation .* complexRegressors;
    regressorsI = zeros(numel(local), numel(pniqComplexSupport));
    regressorsQ = zeros(numel(local), numel(pniqComplexSupport));
    for regressorIndex = 1:numel(pniqComplexSupport)
        descriptor = pniqDescriptors(regressorIndex);
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
        features * pniqCoefficientsI, features * pniqCoefficientsQ);
    prediction = conj(blockRotation) .* predictionNormalized;
    storedPNIQFullPredictions(local, :) = prediction(:, storedColumns);
    pniqFullError = pniqFullError + ...
        sum(abs(prediction - fullSignalTarget(local)).^2, 1);
end
pniqFullNMSE = 10*log10(pniqFullError(:)/targetEnergyFull);

%% Package the two supplementary Ridge families
% Parameters and FLOPs are copied from the principal rows because support is fixed.
variantTargets = repelem(targets, lambdaCount, 1);
variantLambdas = repmat(fixedLambdas, numel(targets), 1);
complexActual = zeros(numel(targets), 1);
complexFLOPs = zeros(numel(targets), 1);
pniqActual = zeros(numel(targets), 1);
pniqFLOPs = zeros(numel(targets), 1);
for targetIndex = 1:numel(targets)
    complexRow = linear.complexTable.TargetRealParameters == targets(targetIndex);
    pniqRow = linear.pniqTable.TargetRealParameters == targets(targetIndex);
    complexActual(targetIndex) = ...
        linear.complexTable.ActualRealParameters(complexRow);
    complexFLOPs(targetIndex) = linear.complexTable.FLOPsPerSample(complexRow);
    pniqActual(targetIndex) = ...
        linear.pniqTable.ActualRealParameters(pniqRow);
    pniqFLOPs(targetIndex) = linear.pniqTable.FLOPsPerSample(pniqRow);
end
complexActual = repelem(complexActual, lambdaCount, 1);
complexFLOPs = repelem(complexFLOPs, lambdaCount, 1);
pniqActual = repelem(pniqActual, lambdaCount, 1);
pniqFLOPs = repelem(pniqFLOPs, lambdaCount, 1);

Model = [repmat(cfg.names.complexGMPDOMP, variantCount, 1); ...
    repmat(cfg.names.pniqGMP, variantCount, 1)];
TargetRealParameters = [variantTargets; variantTargets];
ActualRealParameters = [complexActual; pniqActual];
FixedLambda = [variantLambdas; variantLambdas];
IdentificationNMSEdB = [complexIdentificationNMSE; pniqIdentificationNMSE];
FullSignalNMSEdB = [complexFullNMSE; pniqFullNMSE];
FLOPsPerSample = [complexFLOPs; pniqFLOPs];
MaxAbsRealParameter = [complexMaxAbs; pniqMaxAbs];
fixed.table = table(Model, TargetRealParameters, ActualRealParameters, ...
    FixedLambda, IdentificationNMSEdB, FullSignalNMSEdB, FLOPsPerSample, ...
    MaxAbsRealParameter);
fixed.table.Properties.UserData = struct( ...
    'coefficientRangeDefinition', ...
        string(cfg.sweep.coefficientRangeDefinition), ...
    'linearIdentificationScope', ...
        string(cfg.sweep.linearIdentificationScope), ...
    'linearPrincipalLambda', double(cfg.sweep.linearPrincipalLambda), ...
    'linearLambdaSelection', string(cfg.sweep.linearLambdaSelection), ...
    'fixedRidgeSupportPolicy', ...
        string(cfg.sweep.fixedRidgeSupportPolicy));
fixed.paths = linear.paths;
fixed.pniqPathMap = linear.pniqPathMap;
if storePredictions
    fixed.predictions = struct( ...
        'complexFull', storedComplexFullPredictions, ...
        'pniqFull', storedPNIQFullPredictions);
end
end
