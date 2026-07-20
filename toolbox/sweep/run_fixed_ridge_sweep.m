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
inputRMS = sqrt(mean(abs(x(identificationRows)).^2));
outputRMS = sqrt(mean(abs(y(identificationRows)).^2));

%% Complex GMP: normalize, fit all fixed lambdas, and predict
% The complex column norm is computed once on identification and undone in h.
fprintf('[Fixed ridge] Building one Complex GMP identification matrix...\n');
complexSupport = complexPath(1:maximumFeatures);
identificationU = buildGMPRegressorRows( ...
    x, identificationRows, manager, complexSupport);
complexColumnNorms = sqrt(sum(abs(identificationU).^2, 1)).';
complexCoefficientScales = zeros(maximumFeatures, 1);
for featureIndex = 1:maximumFeatures
    regressor = manager.regPopulation(complexSupport(featureIndex));
    degree = numel(regressor.X) + numel(regressor.Xconj) + ...
        numel(regressor.Xenv);
    complexCoefficientScales(featureIndex) = inputRMS^degree/outputRMS;
end
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
complexMaxAbs = zeros(variantCount, 1);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    for lambdaIndex = 1:lambdaCount
        values = complexCoefficients(1:count, columns(lambdaIndex)) .* ...
            complexCoefficientScales(1:count);
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
fprintf('[Fixed ridge] Building one PN-IQ identification matrix...\n');
selectedMetadata = linear.pnPathMap(1:maximumFeatures, :);
pnComplexSupport = unique(selectedMetadata.SourceRegressorIndex, 'stable');
pnCoefficientScales = zeros(maximumFeatures, 1);
for featureIndex = 1:maximumFeatures
    sourceIndex = selectedMetadata.SourceRegressorIndex(featureIndex);
    regressor = manager.regPopulation(sourceIndex);
    degree = numel(regressor.X) + numel(regressor.Xconj) + ...
        numel(regressor.Xenv);
    pnCoefficientScales(featureIndex) = inputRMS^degree/outputRMS;
end
[~, selectedColumns] = ismember( ...
    selectedMetadata.SourceRegressorIndex, pnComplexSupport);
isQ = selectedMetadata.IsQ;
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
pnMaxAbs = zeros(variantCount, 1);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    for lambdaIndex = 1:lambdaCount
        column = columns(lambdaIndex);
        equivalentCoefficientsI = pnCoefficientsI(1:count, column) .* ...
            pnCoefficientScales(1:count);
        equivalentCoefficientsQ = pnCoefficientsQ(1:count, column) .* ...
            pnCoefficientScales(1:count);
        pnMaxAbs(column) = max(abs([ ...
            equivalentCoefficientsI; equivalentCoefficientsQ]));
    end
end

pnFullError = zeros(1, variantCount);
storedPNFullPredictions = ...
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
MaxAbsRealParameter = [complexMaxAbs; pnMaxAbs];
fixed.table = table(Model, TargetRealParameters, ActualRealParameters, ...
    FixedLambda, IdentificationNMSEdB, FullSignalNMSEdB, FLOPsPerSample, ...
    MaxAbsRealParameter);
fixed.table.Properties.UserData = struct( ...
    'coefficientRangeDefinition', ...
    string(cfg.sweep.coefficientRangeDefinition));
if storePredictions
    fixed.predictions = struct( ...
        'complexFull', storedComplexFullPredictions, ...
        'pnFull', storedPNFullPredictions);
end
end
