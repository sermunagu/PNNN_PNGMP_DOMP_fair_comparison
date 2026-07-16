function selected = buildSelectedFixedLambdaPredictions( ...
    x, y, split, cfg, linear, fixed, selectedParameters)
% buildSelectedFixedLambdaPredictions - Refit six stored-support ridge models.
% One selected target reuses the exact Complex GMP and PN-DOMP support prefixes.
% The full 294-point sweep, support selection, and PNNN are not executed.

requiredLinear = {'complexTable','pnTable','supports','paths','featureMetadata'};
if ~isstruct(linear) || ~all(isfield(linear, requiredLinear)) || ...
        ~isstruct(fixed) || ~isfield(fixed, 'table') || ...
        ~istable(fixed.table)
    error('buildSelectedFixedLambdaPredictions:InvalidArtifacts', ...
        'Compatible linear and fixed-lambda artifacts are required.');
end
complexIndex = uniqueTargetIndex(linear.complexTable, selectedParameters);
pnIndex = uniqueTargetIndex(linear.pnTable, selectedParameters);

selectedLinear = linear;
selectedLinear.complexTable = linear.complexTable(complexIndex, :);
selectedLinear.pnTable = linear.pnTable(pnIndex, :);
selectedLinear.supports.complex = linear.supports.complex(complexIndex);
selectedLinear.supports.pnFeatures = linear.supports.pnFeatures(pnIndex);
selectedLinear.supports.pnComplex = linear.supports.pnComplex(pnIndex);
selectedCfg = cfg;
selectedCfg.sweep.parameterGrid = selectedParameters;

timer = tic;
evaluated = runFixedLambdaLinearSweep( ...
    x, y, split, selectedCfg, selectedLinear, true);
elapsedSeconds = toc(timer);
expected = orderedRows(fixed.table, selectedParameters);
actual = orderedRows(evaluated.table, selectedParameters);
if any(actual.ActualRealParameters ~= expected.ActualRealParameters) || ...
        any(actual.FLOPsPerSample ~= expected.FLOPsPerSample)
    error('buildSelectedFixedLambdaPredictions:ComplexityMismatch', ...
        'Selected ridge complexity no longer matches the saved sweep.');
end
identificationDifference = abs( ...
    actual.IdentificationNMSEdB - expected.IdentificationNMSEdB);
fullSignalDifference = abs( ...
    actual.FullSignalNMSEdB - expected.FullSignalNMSEdB);
if any(identificationDifference > 1e-9) || ...
        any(fullSignalDifference > 1e-9)
    error('buildSelectedFixedLambdaPredictions:NMSEMismatch', ...
        'Selected ridge predictions do not reproduce the saved sweep NMSE.');
end

selected.fixedLambdaComparisonTable = expected;
selected.fullSignalPredictions = struct( ...
    'complexGMP', predictionFamily(evaluated.predictions.complexFull), ...
    'pnIQ', predictionFamily(evaluated.predictions.pnFull));
selected.supports = evaluated.supports;
selected.elapsedSeconds = elapsedSeconds;
selected.metadata = struct('dompInvocationCount', 0, ...
    'pnnnTrainingCount', 0, ...
    'maximumIdentificationNMSEDifferenceDb', ...
        max(identificationDifference), ...
    'maximumFullSignalNMSEDifferenceDb', max(fullSignalDifference));
end

function index = uniqueTargetIndex(value, target)
rows = value.TargetRealParameters == target;
if nnz(rows) ~= 1
    error('buildSelectedFixedLambdaPredictions:MissingTarget', ...
        'The selected target must occur exactly once per linear family.');
end
index = find(rows);
end

function rows = orderedRows(value, target)
models = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"];
lambdas = [1e-3; 1e-4; 1e-5];
rows = value([],:);
for model = models.'
    for lambda = lambdas.'
        match = string(value.Model) == model & ...
            value.TargetRealParameters == target & ...
            value.FixedLambda == lambda;
        if nnz(match) ~= 1
            error('buildSelectedFixedLambdaPredictions:MissingFixedRow', ...
                'Each selected family/lambda row must occur exactly once.');
        end
        rows = [rows; value(match,:)]; %#ok<AGROW>
    end
end
end

function predictions = predictionFamily(values)
if size(values, 2) ~= 3
    error('buildSelectedFixedLambdaPredictions:PredictionCount', ...
        'Each selected family must contain exactly three ridge predictions.');
end
predictions = struct('lambda1e3', values(:, 1), ...
    'lambda1e4', values(:, 2), 'lambda1e5', values(:, 3));
end
