function model = fit_independent_pniq_domp( ...
    x, y, split, cfg, manager, population)
% fit_independent_pniq_domp - Fit the complete independent PN-IQ sweep.
% The file shows phase rotation, I/Q features, PN-DOMP, independent fits,
% phase restoration, predictions, metrics, parameters, and FLOPs.

targets = double(cfg.sweep.parameterGrid(:).');
featureCounts = targets/2;
maximumFeatures = max(featureCounts);
trainRows = split.internalTrainIndices(:);
validationRows = split.internalValidationIndices(:);
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
Component = [repmat("I", nRegressors, 1); ...
    repmat("Q", nRegressors, 1)];
Signature = strings(2*nRegressors, 1);
StructuralZero = false(2*nRegressors, 1);
CanonicalGMP = false(2*nRegressors, 1);
ExactAuxiliaryFallback = false(2*nRegressors, 1);
for index = 1:nRegressors
    descriptor = descriptors(index);
    Signature(index) = descriptor.iSignature;
    Signature(nRegressors + index) = descriptor.qSignature;
    StructuralZero(index) = descriptor.IColumnStructurallyZero;
    StructuralZero(nRegressors + index) = ...
        descriptor.QColumnStructurallyZero;
    CanonicalGMP([index, nRegressors + index]) = descriptor.canonicalGMP;
    ExactAuxiliaryFallback([index, nRegressors + index]) = ...
        ~descriptor.canonicalGMP;
end
RelationSign = ones(2*nRegressors, 1);
rawMetadata = table(SourceRegressorIndex, Component, Signature, ...
    StructuralZero, RelationSign, CanonicalGMP, ExactAuxiliaryFallback);
keptFeatures = find(~StructuralZero);
featureMetadata = rawMetadata(keptFeatures, :);
rawFeatureCount = height(rawMetadata);
effectiveFeatureCount = height(featureMetadata);

%% Internal PN-DOMP path and lambda selection
trainInput = x(trainRows);
trainRotation = complex(ones(size(trainInput)));
nonzero = abs(trainInput) ~= 0;
trainRotation(nonzero) = ...
    conj(trainInput(nonzero)) ./ abs(trainInput(nonzero));
trainTarget = trainRotation .* y(trainRows);

validationInput = x(validationRows);
validationRotation = complex(ones(size(validationInput)));
nonzero = abs(validationInput) ~= 0;
validationRotation(nonzero) = ...
    conj(validationInput(nonzero)) ./ abs(validationInput(nonzero));
validationTarget = y(validationRows);

fprintf('[Linear] Building PN-IQ internal matrices...\n');
trainFeatures = zeros(numel(trainRows), effectiveFeatureCount);
for first = 1:cfg.sweep.candidateBlockSize:numel(trainRows)
    local = first:min(first + cfg.sweep.candidateBlockSize - 1, ...
        numel(trainRows));
    raw = buildFeatures(x, trainRows(local), trainRotation(local), ...
        manager, population, descriptors);
    trainFeatures(local, :) = raw(:, keptFeatures);
end
validationFeatures = zeros(numel(validationRows), effectiveFeatureCount);
for first = 1:cfg.sweep.candidateBlockSize:numel(validationRows)
    local = first:min(first + cfg.sweep.candidateBlockSize - 1, ...
        numel(validationRows));
    raw = buildFeatures(x, validationRows(local), ...
        validationRotation(local), manager, population, descriptors);
    validationFeatures(local, :) = raw(:, keptFeatures);
end

fprintf('[Linear] Computing one PN-IQ PN-DOMP path on internal train...\n');
[trainPath, ~] = selectDOMPSupport( ...
    trainFeatures, trainTarget, maximumFeatures, cfg.gmp.dompOptions);
trainPath = trainPath(:);

selectedLambdas = zeros(size(featureCounts));
validationNMSE = zeros(size(featureCounts));
for targetIndex = 1:numel(targets)
    support = trainPath(1:featureCounts(targetIndex));
    selectedTrain = trainFeatures(:, support);
    featureNorms = sqrt(sum(selectedTrain.^2, 1)).';
    normalizedFeatures = selectedTrain ./ featureNorms.';
    gram = normalizedFeatures.' * normalizedFeatures;
    rhsI = normalizedFeatures.' * real(trainTarget);
    rhsQ = normalizedFeatures.' * imag(trainTarget);
    rankTolerance = max(size(normalizedFeatures)) * ...
        eps(norm(normalizedFeatures, 2));
    candidateNMSE = zeros(numel(cfg.lambdaGrid), 1);

    for lambdaIndex = 1:numel(cfg.lambdaGrid)
        lambda = cfg.lambdaGrid(lambdaIndex);
        if lambda == 0
            normalizedI = lsqminnorm( ...
                normalizedFeatures, real(trainTarget), rankTolerance);
            normalizedQ = lsqminnorm( ...
                normalizedFeatures, imag(trainTarget), rankTolerance);
        else
            regularizedGram = gram + lambda*eye(numel(support));
            normalizedI = regularizedGram \ rhsI;
            normalizedQ = regularizedGram \ rhsQ;
        end
        coefficientsI = normalizedI ./ featureNorms;
        coefficientsQ = normalizedQ ./ featureNorms;
        rotatedPrediction = complex( ...
            validationFeatures(:, support) * coefficientsI, ...
            validationFeatures(:, support) * coefficientsQ);
        prediction = conj(validationRotation) .* rotatedPrediction;
        candidateNMSE(lambdaIndex) = ...
            nmseComplexDb(validationTarget, prediction);
    end

    [validationNMSE(targetIndex), selected] = min(candidateNMSE);
    selectedLambdas(targetIndex) = cfg.lambdaGrid(selected);
end
clear trainFeatures validationFeatures

%% Final PN-DOMP path and independent I/Q fits
identificationInput = x(identificationRows);
identificationRotation = complex(ones(size(identificationInput)));
nonzero = abs(identificationInput) ~= 0;
identificationRotation(nonzero) = ...
    conj(identificationInput(nonzero)) ./ abs(identificationInput(nonzero));
identificationTarget = y(identificationRows);
rotatedIdentificationTarget = ...
    identificationRotation .* identificationTarget;
fullSignalTarget = y(fullSignalRows);

fprintf('[Linear] Building the PN-IQ identification matrix...\n');
identificationFeatures = zeros( ...
    numel(identificationRows), effectiveFeatureCount);
for first = 1:cfg.sweep.candidateBlockSize:numel(identificationRows)
    local = first:min(first + cfg.sweep.candidateBlockSize - 1, ...
        numel(identificationRows));
    raw = buildFeatures(x, identificationRows(local), ...
        identificationRotation(local), manager, population, descriptors);
    identificationFeatures(local, :) = raw(:, keptFeatures);
end

fprintf('[Linear] Computing one PN-IQ PN-DOMP path on identification...\n');
[identificationPath, ~] = selectDOMPSupport( ...
    identificationFeatures, rotatedIdentificationTarget, ...
    maximumFeatures, cfg.gmp.dompOptions);
identificationPath = identificationPath(:);

coefficientsI = zeros(maximumFeatures, numel(targets));
coefficientsQ = zeros(maximumFeatures, numel(targets));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    support = identificationPath(1:count);
    selectedFeatures = identificationFeatures(:, support);
    featureNorms = sqrt(sum(selectedFeatures.^2, 1)).';
    normalizedFeatures = selectedFeatures ./ featureNorms.';
    gram = normalizedFeatures.' * normalizedFeatures;
    rhsI = normalizedFeatures.' * real(rotatedIdentificationTarget);
    rhsQ = normalizedFeatures.' * imag(rotatedIdentificationTarget);
    rankTolerance = max(size(normalizedFeatures)) * ...
        eps(norm(normalizedFeatures, 2));
    lambda = selectedLambdas(targetIndex);

    if lambda == 0
        normalizedI = lsqminnorm(normalizedFeatures, ...
            real(rotatedIdentificationTarget), rankTolerance);
        normalizedQ = lsqminnorm(normalizedFeatures, ...
            imag(rotatedIdentificationTarget), rankTolerance);
    else
        regularizedGram = gram + lambda*eye(count);
        normalizedI = regularizedGram \ rhsI;
        normalizedQ = regularizedGram \ rhsQ;
    end
    coefficientsI(1:count, targetIndex) = normalizedI ./ featureNorms;
    coefficientsQ(1:count, targetIndex) = normalizedQ ./ featureNorms;
end

selectedIdentificationFeatures = identificationFeatures(:, ...
    identificationPath(1:maximumFeatures));
rotatedIdentificationPrediction = complex( ...
    selectedIdentificationFeatures * coefficientsI, ...
    selectedIdentificationFeatures * coefficientsQ);
identificationPredictions = ...
    conj(identificationRotation) .* rotatedIdentificationPrediction;

%% Full-signal prediction and phase restoration
selectedMetadata = featureMetadata( ...
    identificationPath(1:maximumFeatures), :);
complexSupport = unique(selectedMetadata.SourceRegressorIndex, 'stable');
[~, selectedColumns] = ismember( ...
    selectedMetadata.SourceRegressorIndex, complexSupport);
isQ = string(selectedMetadata.Component) == "Q";
selectedColumns(isQ) = selectedColumns(isQ) + numel(complexSupport);

fullPredictions = complex(zeros(numel(fullSignalRows), numel(targets)));
fullBuildCount = 0;
for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, ...
        numel(fullSignalRows));
    xRows = x(fullSignalRows(local));
    rotation = complex(ones(size(xRows)));
    nonzero = abs(xRows) ~= 0;
    rotation(nonzero) = conj(xRows(nonzero)) ./ abs(xRows(nonzero));

    raw = buildFeatures(x, fullSignalRows(local), rotation, ...
        manager, complexSupport, descriptors);
    features = raw(:, selectedColumns);
    rotatedPrediction = complex( ...
        features * coefficientsI, features * coefficientsQ);
    fullPredictions(local, :) = conj(rotation) .* rotatedPrediction;
    fullBuildCount = fullBuildCount + 1;
end

%% Supports, NMSE, parameters, and FLOPs
schema = cfg.sweep.resultSchema;
resultTable = table('Size', [0 numel(schema.names)], ...
    'VariableTypes', schema.types, 'VariableNames', schema.names);
featureSupports = cell(numel(targets), 1);
complexSupports = cell(numel(targets), 1);
for targetIndex = 1:numel(targets)
    features = identificationPath(1:featureCounts(targetIndex));
    metadata = featureMetadata(features, :);
    support = unique(metadata.SourceRegressorIndex, 'stable');
    featureSupports{targetIndex} = features;
    complexSupports{targetIndex} = support;
    reduction = struct('effectiveFeatureCount', featureCounts(targetIndex));
    operations = countModelOperations( ...
        manager.regPopulation, support, reduction);
    cost = countModelFLOPs(operations(4, :), getFLOPConvention());
    values = {"Independent PN-IQ PN-DOMP sweep", "Sweep point", ...
        targets(targetIndex), double(cost.NumRealParameters), ...
        "Validation-selected Ridge", selectedLambdas(targetIndex), ...
        validationNMSE(targetIndex), ...
        nmseComplexDb(identificationTarget, ...
        identificationPredictions(:, targetIndex)), ...
        nmseComplexDb(fullSignalTarget, fullPredictions(:, targetIndex)), ...
        double(cost.FLOPsPerSample), targets(targetIndex), 0, NaN, ...
        "Not applicable", rawFeatureCount, effectiveFeatureCount, ...
        featureCounts(targetIndex), NaN, "Not applicable", NaN, NaN, ...
        NaN, "linear_sweep.mat"};
    resultTable(targetIndex, :) = ...
        cell2table(values, 'VariableNames', schema.names);
end

model.table = resultTable;
model.supports = struct('features', {featureSupports}, ...
    'complex', {complexSupports});
model.paths = struct('train', trainPath, ...
    'identification', identificationPath);
model.lambdas = selectedLambdas;
model.predictions = struct('identification', identificationPredictions, ...
    'full', fullPredictions);
model.featureMetadata = featureMetadata;
model.metadata = struct('candidateFeatures', rawFeatureCount, ...
    'retainedFeatures', effectiveFeatureCount, 'dompInternalTrain', 1, ...
    'dompIdentification', 1, 'matrixInternalTrain', 1, ...
    'matrixInternalValidation', 1, 'matrixIdentification', 1, ...
    'matrixFullSignal', 1, 'fullSignalBuildCount', fullBuildCount);
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
            envelope = envelope .* ...
                abs(x(envelopeRows)).^descriptor.envelopePowers(termIndex);
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
