function model = fit_independent_pniq_domp(x, y, split, cfg, manager, population)
% fit_independent_pniq_domp - Fit the complete independent PN-IQ sweep.
% The file shows phase rotation, I/Q features, PN-DOMP, independent fits,
% phase restoration, predictions, metrics, parameters, and FLOPs.

targets = double(cfg.sweep.parameterGrid(:).');
featureCounts = targets/2;
maximumFeatures = max(featureCounts);
identificationRows = split.identificationIndices(:);
fullSignalRows = split.fullSignalIndices(:);

%% Describe the candidate PN-IQ population
nRegressors = numel(population);
descriptors = repmat(factorizeGMPRegressor( ...
    manager.regPopulation(population(1)), population(1)), nRegressors, 1);

for index = 1:nRegressors
    descriptors(index) = factorizeGMPRegressor( ...
        manager.regPopulation(population(index)), population(index));
end

SourceRegressorIndex = [population; population];
IsQ = [false(nRegressors, 1); true(nRegressors, 1)];
structuralZero = false(2*nRegressors, 1);

for index = 1:nRegressors
    descriptor = descriptors(index);
    structuralZero(nRegressors + index) = descriptor.QColumnStructurallyZero;
end

keptFeatures = find(~structuralZero);
pnFeatureMap = table(SourceRegressorIndex(keptFeatures), ...
    IsQ(keptFeatures), ...
    'VariableNames', {'SourceRegressorIndex','IsQ'});
effectiveFeatureCount = height(pnFeatureMap);

%% Identification PN-DOMP path and independent principal I/Q fits
identificationInput = x(identificationRows);
identificationRotation = complex(ones(size(identificationInput)));
nonzero = abs(identificationInput) ~= 0;
identificationRotation(nonzero) = ...
    conj(identificationInput(nonzero)) ./ abs(identificationInput(nonzero));
identificationTarget = y(identificationRows);
rotatedIdentificationTarget = identificationRotation .* identificationTarget;
fullSignalTarget = y(fullSignalRows);

fprintf('[Linear] Building the PN-IQ identification matrix...\n');
identificationFeatures = zeros(numel(identificationRows), effectiveFeatureCount);

for first = 1:cfg.sweep.candidateBlockSize:numel(identificationRows)
    local = first:min(first + cfg.sweep.candidateBlockSize - 1, numel(identificationRows));
    raw = buildFeatures(x, identificationRows(local), ...
        identificationRotation(local), manager, population, descriptors);
    identificationFeatures(local, :) = raw(:, keptFeatures);
end

fprintf('[Linear] Computing one PN-IQ PN-DOMP path on identification...\n');
identificationPath = selectDOMPSupport( ...
    identificationFeatures, rotatedIdentificationTarget, ...
    maximumFeatures, cfg.gmp.dompOptions.columnTolerance);
identificationPath = identificationPath(:);

% Unit-peak input gives the same unit-norm features because its global scale
% cancels when each homogeneous GMP feature is normalized.
outputPeak = max(abs(y(identificationRows)));

coefficientsI = zeros(maximumFeatures, numel(targets));
coefficientsQ = zeros(maximumFeatures, numel(targets));
comparisonI = coefficientsI;
comparisonQ = coefficientsQ;

for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    support = identificationPath(1:count);
    selectedFeatures = identificationFeatures(:, support);
    featureNorms = sqrt(sum(selectedFeatures.^2, 1)).';
    normalizedFeatures = selectedFeatures ./ featureNorms.';
    rankTolerance = max(size(normalizedFeatures)) * ...
        eps(norm(normalizedFeatures, 2));
    normalizedI = lsqminnorm(normalizedFeatures, ...
        real(rotatedIdentificationTarget), rankTolerance);
    normalizedQ = lsqminnorm(normalizedFeatures, ...
        imag(rotatedIdentificationTarget), rankTolerance);
    comparisonI(1:count, targetIndex) = normalizedI / outputPeak;
    comparisonQ(1:count, targetIndex) = normalizedQ / outputPeak;
    coefficientsI(1:count, targetIndex) = normalizedI ./ featureNorms;
    coefficientsQ(1:count, targetIndex) = normalizedQ ./ featureNorms;
end

selectedIdentificationFeatures = identificationFeatures(:, identificationPath(1:maximumFeatures));
rotatedIdentificationPrediction = complex( ...
    selectedIdentificationFeatures * coefficientsI, ...
    selectedIdentificationFeatures * coefficientsQ);
identificationPredictions = conj(identificationRotation) .* rotatedIdentificationPrediction;

%% Full-signal prediction and phase restoration
selectedMetadata = pnFeatureMap(identificationPath(1:maximumFeatures), :);

complexSupport = unique(selectedMetadata.SourceRegressorIndex, 'stable');

[~, selectedColumns] = ismember(selectedMetadata.SourceRegressorIndex, complexSupport);
isQ = selectedMetadata.IsQ;
selectedColumns(isQ) = selectedColumns(isQ) + numel(complexSupport);

fullPredictions = complex(zeros(numel(fullSignalRows), numel(targets)));

for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    xRows = x(fullSignalRows(local));
    rotation = complex(ones(size(xRows)));
    nonzero = abs(xRows) ~= 0;
    rotation(nonzero) = conj(xRows(nonzero)) ./ abs(xRows(nonzero));

    raw = buildFeatures(x, fullSignalRows(local), rotation, ...
        manager, complexSupport, descriptors);
    features = raw(:, selectedColumns);
    rotatedPrediction = complex(features * coefficientsI, features * coefficientsQ);
    fullPredictions(local, :) = conj(rotation) .* rotatedPrediction;
end

%% NMSE, parameters, and FLOPs
Model = repmat("Independent PN-IQ PN-DOMP sweep", numel(targets), 1);
TargetRealParameters = targets(:);
ActualRealParameters = zeros(numel(targets), 1);
SelectedLambda = zeros(numel(targets), 1);
InternalValidationNMSEdB = nan(numel(targets), 1);
IdentificationNMSEdB = zeros(numel(targets), 1);
FullSignalNMSEdB = zeros(numel(targets), 1);
FLOPsPerSample = zeros(numel(targets), 1);
ActiveWeights = targets(:);
ActiveBiases = zeros(numel(targets), 1);
WeightSparsityPercent = nan(numel(targets), 1);
FineTuneEpochs = nan(numel(targets), 1);
MaxAbsRealParameter = zeros(numel(targets), 1);

for targetIndex = 1:numel(targets)
    features = identificationPath(1:featureCounts(targetIndex));
    metadata = pnFeatureMap(features, :);
    support = unique(metadata.SourceRegressorIndex, 'stable');
    operations = countModelOperations( ...
        manager.regPopulation, support, featureCounts(targetIndex));
    cost = countModelFLOPs(operations(4, :), getFLOPConvention());
    ActualRealParameters(targetIndex) = double(cost.NumRealParameters);
    IdentificationNMSEdB(targetIndex) = nmseComplexDb( ...
        identificationTarget, identificationPredictions(:, targetIndex));
    FullSignalNMSEdB(targetIndex) = nmseComplexDb( ...
        fullSignalTarget, fullPredictions(:, targetIndex));
    FLOPsPerSample(targetIndex) = double(cost.FLOPsPerSample);
    count = featureCounts(targetIndex);
    MaxAbsRealParameter(targetIndex) = max(abs([ ...
        comparisonI(1:count, targetIndex); ...
        comparisonQ(1:count, targetIndex)]));
end

resultTable = table(Model, TargetRealParameters, ActualRealParameters, ...
    SelectedLambda, InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, ActiveWeights, ActiveBiases, ...
    WeightSparsityPercent, FineTuneEpochs, MaxAbsRealParameter);

model.table = resultTable;
model.path = identificationPath;
model.fullPredictions = fullPredictions;
model.pnPathMap = pnFeatureMap(identificationPath, :);
end






function features = buildFeatures( ...
    x, rows, rotation, manager, support, descriptors)
% Build the real I/Q representation of phase-normalized GMP regressors.

complexRegressors = buildGMPRegressorRows(x, rows, manager, support);
phaseNormalized = rotation .* complexRegressors;
regressorsI = zeros(numel(rows), numel(support));
regressorsQ = zeros(numel(rows), numel(support));
signalLength = numel(x);

for localIndex = 1:numel(support)
    descriptor = descriptors(support(localIndex));
    if descriptor.canonicalGMP
        carrierRows = mod( ...
            rows - descriptor.carrierLag - 1, signalLength) + 1;
        normalizedCarrier = rotation .* x(carrierRows);
        envelope = ones(numel(rows), 1);
        for termIndex = 1:numel(descriptor.envelopeLags)
            envelopeRows = mod(rows - ...
                descriptor.envelopeLags(termIndex) - 1, signalLength) + 1;
            envelope = envelope .* abs(x(envelopeRows)).^ ...
                descriptor.envelopePowers(termIndex);
        end
        regressorsI(:, localIndex) = real(normalizedCarrier) .* envelope;
        regressorsQ(:, localIndex) = imag(normalizedCarrier) .* envelope;
        if descriptor.QColumnStructurallyZero
            regressorsQ(:, localIndex) = 0;
        end
    else
        regressorsI(:, localIndex) = real(phaseNormalized(:, localIndex));
        regressorsQ(:, localIndex) = imag(phaseNormalized(:, localIndex));
    end
end
features = [regressorsI, regressorsQ];
end
