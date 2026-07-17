function model = fit_complex_gmp_domp(x, y, split, cfg, manager, population)
% fit_complex_gmp_domp - Fit the complete Complex GMP-DOMP sweep.
% DOMP and lambda selection use the internal split; final coefficients use
% identification only, and the stored predictions cover identification/full.

targets = double(cfg.sweep.parameterGrid(:).');
featureCounts = targets/2;
maximumFeatures = max(featureCounts);
trainRows = split.internalTrainIndices(:);
validationRows = split.internalValidationIndices(:);
identificationRows = split.identificationIndices(:);
fullSignalRows = split.fullSignalIndices(:);

%% Internal DOMP path and lambda selection
fprintf('[Linear] Building Complex GMP internal matrices...\n');

trainU = buildGMPRegressorRows(x, trainRows, manager, population);
validationU = buildGMPRegressorRows(x, validationRows, manager, population);

trainTarget = y(trainRows);
validationTarget = y(validationRows);

fprintf('[Linear] Computing one Complex GMP DOMP path on internal train...\n');

[trainPath, ~] = selectDOMPSupport(trainU, trainTarget, maximumFeatures, cfg.gmp.dompOptions);
trainPath = trainPath(:);
trainNorms = sqrt(sum(abs(trainU).^2, 1)).';

selectedLambdas = zeros(size(featureCounts));
validationNMSE = zeros(size(featureCounts));

for targetIndex = 1:numel(targets)
    support = trainPath(1:featureCounts(targetIndex));
    columnNorms = trainNorms(support);
    normalizedU = trainU(:, support) ./ columnNorms.';
    gram = normalizedU' * normalizedU;
    rhs = normalizedU' * trainTarget;
    rankTolerance = max(size(normalizedU))*eps(norm(normalizedU, 2));
    candidateNMSE = zeros(numel(cfg.lambdaGrid), 1);

    for lambdaIndex = 1:numel(cfg.lambdaGrid)
        lambda = cfg.lambdaGrid(lambdaIndex);
        if lambda == 0
            normalizedCoefficients = lsqminnorm(normalizedU, trainTarget, rankTolerance);
        else
            normalizedCoefficients = (gram + lambda*eye(numel(support))) \ rhs;
        end
        coefficients = normalizedCoefficients ./ columnNorms;
        prediction = validationU(:, support) * coefficients;
        candidateNMSE(lambdaIndex) = nmseComplexDb(validationTarget, prediction);
    end

    [validationNMSE(targetIndex), selected] = min(candidateNMSE);
    selectedLambdas(targetIndex) = cfg.lambdaGrid(selected);
end
clear trainU validationU

%% Final DOMP path and identification fits
identificationTarget = y(identificationRows);
fullSignalTarget = y(fullSignalRows);
fprintf('[Linear] Building the Complex GMP identification matrix...\n');
identificationU = buildGMPRegressorRows(x, identificationRows, manager, population);

fprintf('[Linear] Computing one Complex GMP DOMP path on identification...\n');
[identificationPath, ~] = selectDOMPSupport(identificationU, identificationTarget, maximumFeatures, cfg.gmp.dompOptions);
identificationPath = identificationPath(:);
identificationNorms = sqrt(sum(abs(identificationU).^2, 1)).';

coefficients = complex(zeros(maximumFeatures, numel(targets)));
for targetIndex = 1:numel(targets)
    count = featureCounts(targetIndex);
    support = identificationPath(1:count);
    columnNorms = identificationNorms(support);
    normalizedU = identificationU(:, support) ./ columnNorms.';
    gram = normalizedU' * normalizedU;
    rhs = normalizedU' * identificationTarget;
    rankTolerance = max(size(normalizedU))*eps(norm(normalizedU, 2));
    lambda = selectedLambdas(targetIndex);

    if lambda == 0
        normalizedCoefficients = lsqminnorm(normalizedU, identificationTarget, rankTolerance);
    else
        normalizedCoefficients = (gram + lambda*eye(count)) \ rhs;
    end
    
    coefficients(1:count, targetIndex) = normalizedCoefficients ./ columnNorms;
end

identificationPredictions = identificationU(:, ...
    identificationPath(1:maximumFeatures)) * coefficients;

%% Full-signal prediction
fprintf('[Linear] Evaluating %d Complex GMP targets on the full signal...\n', numel(targets));
fullPredictions = complex(zeros(numel(fullSignalRows), numel(targets)));
fullBuildCount = 0;

for first = 1:cfg.gmp.blockSize:numel(fullSignalRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullSignalRows));
    
    U = buildGMPRegressorRows(x, fullSignalRows(local), manager, identificationPath(1:maximumFeatures));
    
    fullPredictions(local, :) = U * coefficients;
    fullBuildCount = fullBuildCount + 1;
end

%% Supports, NMSE, parameters, and FLOPs
schema = cfg.sweep.resultSchema;
resultTable = table('Size', [0 numel(schema.names)], 'VariableTypes', schema.types, 'VariableNames', schema.names);
supports = cell(numel(targets), 1);

for targetIndex = 1:numel(targets)
    support = identificationPath(1:featureCounts(targetIndex));
    supports{targetIndex} = support;
    operations = countModelOperations(manager.regPopulation, support);
    cost = countModelFLOPs(operations(1, :), getFLOPConvention());
    values = {"Complex GMP DOMP sweep", "Sweep point", ...
        targets(targetIndex), double(cost.NumRealParameters), ...
        "Validation-selected Ridge", selectedLambdas(targetIndex), ...
        validationNMSE(targetIndex), ...
        nmseComplexDb(identificationTarget, ...
        identificationPredictions(:, targetIndex)), ...
        nmseComplexDb(fullSignalTarget, fullPredictions(:, targetIndex)), ...
        double(cost.FLOPsPerSample), targets(targetIndex), 0, NaN, ...
        "Not applicable", NaN, NaN, featureCounts(targetIndex), NaN, ...
        "Not applicable", NaN, NaN, NaN, "linear_sweep.mat"};
    resultTable(targetIndex, :) = cell2table(values, 'VariableNames', schema.names);
end

model.table = resultTable;
model.supports = supports;
model.paths = struct('train', trainPath, 'identification', identificationPath);
model.lambdas = selectedLambdas;
model.predictions = struct('identification', identificationPredictions, 'full', fullPredictions);
model.metadata = struct('dompInternalTrain', 1, ...
    'dompIdentification', 1, 'matrixInternalTrain', 1, ...
    'matrixInternalValidation', 1, 'matrixIdentification', 1, ...
    'matrixFullSignal', 1, 'fullSignalBuildCount', fullBuildCount);
end
