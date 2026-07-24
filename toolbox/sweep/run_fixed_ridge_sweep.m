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
outputPeak = max(abs(identificationTarget));
if ~isfinite(outputPeak) || outputPeak <= 0
    error('run_fixed_ridge_sweep:InvalidIdentificationOutputPeak', ...
        'The identification output peak must be finite and positive.');
end
unitPeakIdentificationTarget = identificationTarget / outputPeak;

%% Complex GMP: peak-normalize, fit all fixed lambdas, and predict
% Global input scaling cancels when each homogeneous GMP column is peak-normalized.
fprintf('[Fixed ridge] Building one %s identification matrix...\n', ...
    cfg.names.complexGMPDOMP);
complexSupport = complexPath(1:maximumFeatures);
identificationU = buildGMPRegressorRows( ...
    x, identificationRows, manager, complexSupport);
complexPredictionCoefficients = complex(zeros(maximumFeatures, variantCount));
complexUnitPeakRegressionCoefficientPaths = ...
    complex(zeros(maximumFeatures, variantCount));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    selectedRegressors = identificationU(:, 1:count);
    regressorPeaks = max(abs(selectedRegressors), [], 1).';
    invalidColumn = find(~isfinite(regressorPeaks) | ...
        regressorPeaks <= 0, 1);
    if ~isempty(invalidColumn)
        error('run_fixed_ridge_sweep:InvalidComplexRegressorPeak', ...
            ['%s fixed-Ridge target %d selected column %d has a ' ...
            'nonfinite or nonpositive peak.'], cfg.names.complexGMPDOMP, ...
            targets(targetIndex), complexSupport(invalidColumn));
    end
    unitPeakRegressors = selectedRegressors ./ regressorPeaks.';
    gram = unitPeakRegressors' * unitPeakRegressors;
    rhs = unitPeakRegressors' * unitPeakIdentificationTarget;
    for lambdaIndex = 1:lambdaCount
        column = columns(lambdaIndex);
        regularizedGram = gram + ...
            fixedLambdas(lambdaIndex)*eye(count, 'like', gram);
        unitPeakRegressionCoefficients = regularizedGram \ rhs;
        complexUnitPeakRegressionCoefficientPaths(1:count, column) = ...
            unitPeakRegressionCoefficients;
        complexPredictionCoefficients(1:count, column) = ...
            (unitPeakRegressionCoefficients * outputPeak) ./ ...
            regressorPeaks;
    end
end
complexIdentificationPredictions = ...
    identificationU * complexPredictionCoefficients;
complexIdentificationError = sum(abs( ...
    complexIdentificationPredictions - identificationTarget).^2, 1);
complexIdentificationNMSE = ...
    10*log10(complexIdentificationError(:)/targetEnergyIdentification);
complexMaxAbs = zeros(variantCount, 1);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    for lambdaIndex = 1:lambdaCount
        column = columns(lambdaIndex);
        activeUnitPeakRegressionCoefficients = ...
            complexUnitPeakRegressionCoefficientPaths(1:count, column);
        complexMaxAbs(column) = ...
            max(abs(activeUnitPeakRegressionCoefficients(:)));
    end
end

complexFullError = zeros(1, variantCount);
storedComplexFullPredictions = complex(zeros(numel(fullSignalRows), numel(storedColumns)));
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    U = buildGMPRegressorRows(x, fullSignalRows(local), manager, complexSupport);
    prediction = U * complexPredictionCoefficients;
    storedComplexFullPredictions(local, :) = prediction(:, storedColumns);
    complexFullError = complexFullError + sum(abs(prediction - fullSignalTarget(local)).^2, 1);
end
complexFullNMSE = 10*log10(complexFullError(:)/targetEnergyFull);
clear identificationU complexPredictionCoefficients complexIdentificationPredictions

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
phaseNormalizedIdentificationTarget = ...
    identificationRotation .* identificationTarget;
unitPeakPhaseNormalizedIdentificationTarget = ...
    phaseNormalizedIdentificationTarget / outputPeak;

pniqPredictionCoefficientsI = zeros(maximumFeatures, variantCount);
pniqPredictionCoefficientsQ = zeros(maximumFeatures, variantCount);
pniqUnitPeakRegressionCoefficientPathsI = ...
    zeros(maximumFeatures, variantCount);
pniqUnitPeakRegressionCoefficientPathsQ = ...
    zeros(maximumFeatures, variantCount);
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    columns = (targetIndex - 1)*lambdaCount + (1:lambdaCount);
    selectedFeatures = identificationFeatures(:, 1:count);
    featurePeaks = max(abs(selectedFeatures), [], 1).';
    invalidFeature = find(~isfinite(featurePeaks) | featurePeaks <= 0, 1);
    if ~isempty(invalidFeature)
        error('run_fixed_ridge_sweep:InvalidPNIQFeaturePeak', ...
            ['%s fixed-Ridge target %d selected feature %d has a ' ...
            'nonfinite or nonpositive peak.'], cfg.names.pniqGMP, ...
            targets(targetIndex), linear.paths.pniq(invalidFeature));
    end
    unitPeakFeatures = selectedFeatures ./ featurePeaks.';
    gram = unitPeakFeatures.' * unitPeakFeatures;
    rhsI = unitPeakFeatures.' * ...
        real(unitPeakPhaseNormalizedIdentificationTarget);
    rhsQ = unitPeakFeatures.' * ...
        imag(unitPeakPhaseNormalizedIdentificationTarget);
    for lambdaIndex = 1:lambdaCount
        regularizedGram = gram + ...
            fixedLambdas(lambdaIndex)*eye(count, 'like', gram);
        column = columns(lambdaIndex);
        unitPeakRegressionCoefficientsI = ...
            regularizedGram \ rhsI(1:count);
        unitPeakRegressionCoefficientsQ = ...
            regularizedGram \ rhsQ(1:count);
        pniqPredictionCoefficientsI(1:count, column) = ...
            (unitPeakRegressionCoefficientsI * outputPeak) ./ ...
            featurePeaks;
        pniqPredictionCoefficientsQ(1:count, column) = ...
            (unitPeakRegressionCoefficientsQ * outputPeak) ./ ...
            featurePeaks;
        pniqUnitPeakRegressionCoefficientPathsI(1:count, column) = ...
            unitPeakRegressionCoefficientsI;
        pniqUnitPeakRegressionCoefficientPathsQ(1:count, column) = ...
            unitPeakRegressionCoefficientsQ;
    end
end
pniqPhaseNormalizedIdentificationPrediction = ...
    identificationFeatures * pniqPredictionCoefficientsI + ...
    1j * (identificationFeatures * pniqPredictionCoefficientsQ);
pniqIdentificationPredictions = ...
    conj(identificationRotation) .* ...
    pniqPhaseNormalizedIdentificationPrediction;
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
        activeI = pniqUnitPeakRegressionCoefficientPathsI( ...
            1:count, column);
        activeQ = pniqUnitPeakRegressionCoefficientPathsQ( ...
            1:count, column);
        pniqMaxAbs(column) = max(abs([ ...
            activeI(:); activeQ(:)]));
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
    phaseNormalizedPrediction = complex( ...
        features * pniqPredictionCoefficientsI, ...
        features * pniqPredictionCoefficientsQ);
    prediction = conj(blockRotation) .* phaseNormalizedPrediction;
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
