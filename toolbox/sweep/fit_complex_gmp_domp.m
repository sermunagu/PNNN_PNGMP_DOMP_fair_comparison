function model = fit_complex_gmp_domp(x, y, split, cfg, manager, population)
% fit_complex_gmp_domp - Fit the complete Complex GMP-DOMP sweep.
% One DOMP path and all principal least-squares fits use identification.

targets = double(cfg.sweep.parameterGrid(:).');
featureCounts = targets/2;
maximumFeatures = max(featureCounts);
identificationRows = split.identificationIndices(:);
fullSignalRows = split.fullSignalIndices(:);

%% Identification DOMP path and principal least-squares fits
identificationTarget = y(identificationRows);
fullSignalTarget = y(fullSignalRows);
fprintf('[Linear] Building the Complex GMP identification matrix...\n');
identificationU = buildGMPRegressorRows(x, identificationRows, manager, population);

fprintf('[Linear] Computing one Complex GMP DOMP path on identification...\n');
identificationPath = selectDOMPSupport(identificationU, ...
    identificationTarget, maximumFeatures, ...
    cfg.gmp.dompOptions.columnTolerance);
identificationPath = identificationPath(:);
identificationNorms = sqrt(sum(abs(identificationU).^2, 1)).';

% Unit-peak input gives the same unit-norm columns because its global scale
% cancels when each homogeneous GMP column is normalized.
outputPeak = max(abs(y(identificationRows)));

coefficients = complex(zeros(maximumFeatures, numel(targets)));
comparisonCoefficients = coefficients;
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    support = identificationPath(1:count);
    columnNorms = identificationNorms(support);
    normalizedU = identificationU(:, support) ./ columnNorms.';
    rankTolerance = max(size(normalizedU))*eps(norm(normalizedU, 2));
    normalizedCoefficients = lsqminnorm( ...
        normalizedU, identificationTarget, rankTolerance);

    comparisonCoefficients(1:count, targetIndex) = ...
        normalizedCoefficients / outputPeak;
    coefficients(1:count, targetIndex) = normalizedCoefficients ./ columnNorms;
end

identificationPredictions = identificationU(:, ...
    identificationPath(1:maximumFeatures)) * coefficients;

%% Full-signal prediction
fprintf('[Linear] Evaluating %d Complex GMP targets on the full signal...\n', numel(targets));
fullPredictions = complex(zeros(numel(fullSignalRows), numel(targets)));

for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    U = buildGMPRegressorRows(x, fullSignalRows(local), manager, ...
        identificationPath(1:maximumFeatures));
    
    fullPredictions(local, :) = U * coefficients;
end

%% NMSE, parameters, and FLOPs
Model = repmat("Complex GMP DOMP sweep", numel(targets), 1);
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
    support = identificationPath(1:featureCounts(targetIndex));
    operations = countModelOperations(manager.regPopulation, support);
    cost = countModelFLOPs(operations(1, :), getFLOPConvention());
    ActualRealParameters(targetIndex) = double(cost.NumRealParameters);
    IdentificationNMSEdB(targetIndex) = nmseComplexDb( ...
        identificationTarget, identificationPredictions(:, targetIndex));
    FullSignalNMSEdB(targetIndex) = nmseComplexDb( ...
        fullSignalTarget, fullPredictions(:, targetIndex));
    FLOPsPerSample(targetIndex) = double(cost.FLOPsPerSample);
    activeCoefficients = comparisonCoefficients( ...
        1:featureCounts(targetIndex), ...
        targetIndex);
    MaxAbsRealParameter(targetIndex) = max([ ...
        abs(real(activeCoefficients)); ...
        abs(imag(activeCoefficients))]);
end

resultTable = table(Model, TargetRealParameters, ActualRealParameters, ...
    SelectedLambda, InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, ActiveWeights, ActiveBiases, ...
    WeightSparsityPercent, FineTuneEpochs, MaxAbsRealParameter);

model.table = resultTable;
model.path = identificationPath;
model.fullPredictions = fullPredictions;
end
