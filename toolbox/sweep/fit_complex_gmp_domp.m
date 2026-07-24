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
fprintf('[Linear] Building the %s identification matrix...\n', ...
    cfg.names.complexGMPDOMP);
identificationU = buildGMPRegressorRows(x, identificationRows, manager, population);

fprintf(['[Linear] Computing one DOMP support path for %s ' ...
    'on identification...\n'], cfg.names.complexGMPDOMP);
identificationPath = selectDOMPSupport(identificationU, ...
    identificationTarget, maximumFeatures, ...
    cfg.gmp.dompOptions.columnTolerance);
identificationPath = identificationPath(:);

outputPeak = max(abs(identificationTarget));
if ~isfinite(outputPeak) || outputPeak <= 0
    error('fit_complex_gmp_domp:InvalidIdentificationOutputPeak', ...
        'The identification output peak must be finite and positive.');
end
unitPeakIdentificationTarget = identificationTarget / outputPeak;

predictionCoefficients = complex(zeros(maximumFeatures, numel(targets)));
unitPeakRegressionCoefficientPaths = complex(zeros(maximumFeatures, numel(targets)));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    support = identificationPath(1:count);
    selectedRegressors = identificationU(:, support);
    regressorPeaks = max(abs(selectedRegressors), [], 1).';
    invalidColumn = find(~isfinite(regressorPeaks) | ...
        regressorPeaks <= 0, 1);
    if ~isempty(invalidColumn)
        error('fit_complex_gmp_domp:InvalidRegressorPeak', ...
            ['%s target %d selected column %d has a nonfinite or ' ...
            'nonpositive peak.'], cfg.names.complexGMPDOMP, ...
            targets(targetIndex), support(invalidColumn));
    end
    
    unitPeakRegressors = selectedRegressors ./ regressorPeaks.';
    rankTolerance = max(size(unitPeakRegressors)) *  eps(norm(unitPeakRegressors, 2));

    unitPeakRegressionCoefficients = lsqminnorm(unitPeakRegressors, unitPeakIdentificationTarget, rankTolerance);

    unitPeakRegressionCoefficientPaths(1:count, targetIndex) = unitPeakRegressionCoefficients;

    predictionCoefficients(1:count, targetIndex) = ...
        (unitPeakRegressionCoefficients * outputPeak) ./ regressorPeaks;
end

identificationPredictions = identificationU(:, identificationPath(1:maximumFeatures)) * predictionCoefficients;

%% Full-signal prediction
fprintf('[Linear] Evaluating %d %s targets on the full signal...\n', ...
    numel(targets), cfg.names.complexGMPDOMP);
fullPredictions = complex(zeros(numel(fullSignalRows), numel(targets)));

for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    U = buildGMPRegressorRows(x, fullSignalRows(local), manager, ...
        identificationPath(1:maximumFeatures));

    fullPredictions(local, :) = U * predictionCoefficients;
end

%% NMSE, parameters, and FLOPs
Model = repmat(cfg.names.complexGMPDOMP, numel(targets), 1);
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
    IdentificationNMSEdB(targetIndex) = nmseComplexDb(identificationTarget, identificationPredictions(:, targetIndex));
    FullSignalNMSEdB(targetIndex) = nmseComplexDb(fullSignalTarget, fullPredictions(:, targetIndex));
    FLOPsPerSample(targetIndex) = double(cost.FLOPsPerSample);

    activeUnitPeakRegressionCoefficients = ...
        unitPeakRegressionCoefficientPaths( ...
        1:featureCounts(targetIndex), targetIndex);
    MaxAbsRealParameter(targetIndex) = ...
        max(abs(activeUnitPeakRegressionCoefficients(:)));
end

resultTable = table(Model, TargetRealParameters, ActualRealParameters, ...
    SelectedLambda, InternalValidationNMSEdB, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, ActiveWeights, ActiveBiases, ...
    WeightSparsityPercent, FineTuneEpochs, MaxAbsRealParameter);

model.table = resultTable;
model.path = identificationPath;
model.fullPredictions = fullPredictions;
end
