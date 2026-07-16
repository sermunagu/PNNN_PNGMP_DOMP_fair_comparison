function results = run_fair_PNNN_vs_PNGMP_DOMP( ...
    selectedParameters, sweep)
% run_fair_PNNN_vs_PNGMP_DOMP - Compare one signed sweep parameter point.
% Linear and sparse-PNNN predictions are loaded from compatible checkpoints.
% This function performs no support selection, fitting, pruning, or training.

if nargin < 1
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingTarget', ...
        'An explicit active real-parameter target is required.');
end
validateTarget(selectedParameters);
if nargin < 2 || isempty(sweep)
    sweep = run_parameter_sweep(20:10:500);
end

[selectedRows, sweepDirectory] = validateSweep( ...
    sweep, selectedParameters);
signature = sweep.experimentSignature;
identity = sweep.sweepIdentity;
[targetFullSignal, fullSignalIndices, sampleRateHz, sampleRateSource, ...
    fittingContext] = ...
    recoverFullSignalTarget(sweep, signature);

summaryArtifact = loadSignedArtifact(fullfile(sweepDirectory, ...
    'complexity_sweep.mat'), identity, signature);
validateSummary(summaryArtifact.payload, selectedRows, selectedParameters);

linearArtifact = loadSignedArtifact(fullfile(sweepDirectory, ...
    'linear_sweep.mat'), identity, signature);
denseArtifact = loadSignedArtifact(fullfile(sweepDirectory, ...
    'sweep_dense_source.mat'), identity, signature);
pnnnFilename = sprintf('pnnn_target_%04d.mat', selectedParameters);
pnnnArtifact = loadSignedArtifact(fullfile(sweepDirectory, ...
    pnnnFilename), identity, signature);

[linearPredictions, linearSupports] = selectLinearPoint( ...
    linearArtifact.payload, selectedRows, selectedParameters);
pnnnPrediction = selectPNNNPoint(pnnnArtifact.payload, ...
    denseArtifact.payload, selectedRows, selectedParameters);
validatePredictionSizes(linearPredictions, pnnnPrediction, ...
    targetFullSignal, fullSignalIndices);
[fixedLambdaTable, fixedLambdaPredictions, ridgePredictionTime] = ...
    selectedFixedLambdaPredictions(sweep, sweepDirectory, identity, ...
    signature, linearArtifact.payload, fittingContext, selectedParameters);
validateFixedPredictionSizes(fixedLambdaTable, fixedLambdaPredictions, ...
    numel(targetFullSignal));

comparisonTable = selectedRows(:, {'Model','ActualRealParameters', ...
    'FullSignalNMSEdB','FLOPsPerSample'});
comparisonTable.Model = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"; ...
    "Sparse PNNN N12"];

results = struct( ...
    'selectedParameters', selectedParameters, ...
    'comparisonTable', comparisonTable, ...
    'resultDirectory', string(sweepDirectory), ...
    'sweepResultDirectory', string(sweepDirectory), ...
    'reusedLinearSweep', true, ...
    'reusedPNNNPoint', true, ...
    'reusedDenseSource', true, ...
    'selectedSweepRows', selectedRows, ...
    'linearSupports', linearSupports, ...
    'targetFullSignal', targetFullSignal, ...
    'fullSignalIndices', fullSignalIndices, ...
    'sampleRateHz', sampleRateHz, ...
    'sampleRateSource', sampleRateSource, ...
    'fixedLambdaComparisonTable', fixedLambdaTable, ...
    'fixedLambdaFullSignalPredictions', fixedLambdaPredictions, ...
    'ridgePredictionTimeSeconds', ridgePredictionTime, ...
    'fullSignalPredictions', struct( ...
        'complexGMP', linearPredictions.complexGMP, ...
        'pnIQ', linearPredictions.pnIQ, ...
        'sparsePNNNN12', pnnnPrediction));

disp(comparisonTable);
fprintf('Reused linear sweep: YES\n');
fprintf('Reused fixed-lambda checkpoint: YES\n');
fprintf('Reused sparse PNNN target %d: YES\n', selectedParameters);
fprintf('Selected fixed-lambda predictions completed in %.2f s.\n', ...
    ridgePredictionTime);
fprintf('Sweep results: %s\n', sweepDirectory);
end

function [target, indices, sampleRateHz, source, fittingContext] = ...
    recoverFullSignalTarget(sweep, expectedSignature)
provided = {'targetFullSignal','fullSignalIndices','sampleRateHz'};
if all(isfield(sweep, provided))
    target = sweep.targetFullSignal(:);
    indices = double(sweep.fullSignalIndices(:));
    sampleRateHz = double(sweep.sampleRateHz);
    source = "Provided signed-sweep context";
    if isfield(sweep, 'sampleRateSource')
        source = string(sweep.sampleRateSource);
    end
    fittingContext = [];
    validateTargetContext(target, indices, sampleRateHz);
    return;
end

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));
addpath(fullfile(projectRoot, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(projectRoot, 'toolbox', 'splits'));
cfg = getFairDOMPComparisonConfig(projectRoot);
measurement = load(cfg.measurementFile, 'x', 'y', 'fs', 'info_signal');
if ~all(isfield(measurement, {'x','y'}))
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidMeasurement', ...
        'The configured measurement must contain x and y.');
end
[x, y] = selectXYByMapping( ...
    measurement.x, measurement.y, cfg.mappingMode);
x = x(:);
y = y(:);
actualSignature = buildExperimentSignature(x, y, cfg);
if ~isequaln(actualSignature, expectedSignature)
    error('run_fair_PNNN_vs_PNGMP_DOMP:SignatureMismatch', ...
        'The configured measurement does not match the signed sweep.');
end
if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end
split = buildCommonComparisonSplit(x, y, cfg);
indices = double(split.fullSignalIndices(:));
target = y(indices);
[sampleRateHz, source] = measurementSampleRate(measurement);
fittingContext = struct('x', x, 'y', y, 'split', split, 'cfg', cfg);
validateTargetContext(target, indices, sampleRateHz);
end

function [tableValue, predictions, elapsedSeconds] = ...
    selectedFixedLambdaPredictions(sweep, directory, identity, signature, ...
    linear, fittingContext, target)
provided = {'fixedLambdaComparisonTable', ...
    'fixedLambdaFullSignalPredictions'};
if all(isfield(sweep, provided))
    tableValue = sweep.fixedLambdaComparisonTable;
    predictions = sweep.fixedLambdaFullSignalPredictions;
    elapsedSeconds = 0;
    return;
end
if isempty(fittingContext)
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingFittingContext', ...
        'Selected fixed-lambda predictions require the signed signal context.');
end
fixedArtifact = loadSignedArtifact(fullfile(directory, ...
    'fixed_lambda_linear_sweep.mat'), identity, signature);
selected = buildSelectedFixedLambdaPredictions( ...
    fittingContext.x, fittingContext.y, fittingContext.split, ...
    fittingContext.cfg, linear, fixedArtifact.payload, target);
tableValue = selected.fixedLambdaComparisonTable;
predictions = selected.fullSignalPredictions;
elapsedSeconds = selected.elapsedSeconds;
end

function [sampleRateHz, source] = measurementSampleRate(measurement)
if isfield(measurement, 'fs') && isscalar(measurement.fs) && ...
        isfinite(measurement.fs) && measurement.fs > 0
    sampleRateHz = double(measurement.fs);
    source = "measurement.fs";
    if isfield(measurement, 'info_signal') && ...
            isstruct(measurement.info_signal) && ...
            isfield(measurement.info_signal, 'fsovs')
        oversampledRate = double(measurement.info_signal.fsovs);
        if ~isscalar(oversampledRate) || ~isfinite(oversampledRate) || ...
                abs(oversampledRate - sampleRateHz) > ...
                eps(max(oversampledRate, sampleRateHz))
            error('run_fair_PNNN_vs_PNGMP_DOMP:SampleRateMismatch', ...
                'Measurement fs and info_signal.fsovs disagree.');
        end
        source = "measurement.fs, confirmed by info_signal.fsovs";
    end
elseif isfield(measurement, 'info_signal') && ...
        isstruct(measurement.info_signal) && ...
        isfield(measurement.info_signal, 'fsovs')
    sampleRateHz = double(measurement.info_signal.fsovs);
    source = "measurement.info_signal.fsovs";
else
    % This local fallback corresponds to the expected capture rate.
    sampleRateHz = 614.4e6;
    source = "Local fallback for the current capture";
end
end

function validateTargetContext(target, indices, sampleRateHz)
valid = ~isempty(target) && numel(target) == numel(indices) && ...
    all(isfinite(target)) && all(isfinite(indices)) && ...
    all(indices == floor(indices)) && all(indices > 0) && ...
    isscalar(sampleRateHz) && isfinite(sampleRateHz) && sampleRateHz > 0;
if ~valid
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidTargetContext', ...
        'The full-signal target, indices, or sample rate are invalid.');
end
end

function validateTarget(target)
valid = isnumeric(target) && isreal(target) && isscalar(target) && ...
    isfinite(target) && target == fix(target);
if ~valid
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidTarget', ...
        'The parameter target must be a finite real integer scalar.');
end
minimumTarget = 15;
if target < minimumTarget
    error('run_fair_PNNN_vs_PNGMP_DOMP:TargetBelowPNNNMinimum', ...
        ['The target must include 14 protected N12 biases and at least ' ...
        'one active weight.']);
end
end

function [rows, directory] = validateSweep(sweep, target)
required = {'results','resultDirectory','experimentSignature','sweepIdentity'};
if ~isstruct(sweep) || ~all(isfield(sweep, required)) || ...
        ~istable(sweep.results) || ...
        ~validSignature(sweep.experimentSignature) || ...
        ~isstruct(sweep.sweepIdentity) || ...
        ~isfield(sweep.sweepIdentity, 'digest')
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidSweep', ...
        'A signed sweep result structure is required.');
end
directory = normalizePath(sweep.resultDirectory);
if ~isfolder(directory)
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingSweepDirectory', ...
        'The signed sweep result directory does not exist.');
end

requiredColumns = {'Model','SweepRole','TargetRealParameters', ...
    'ActualRealParameters','FullSignalNMSEdB','FLOPsPerSample', ...
    'ActiveWeights','ActiveBiases'};
if ~all(ismember(requiredColumns, sweep.results.Properties.VariableNames))
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidSweepTable', ...
        'The sweep table is missing required comparison columns.');
end
expectedModels = ["Complex GMP DOMP sweep"; ...
    "Independent PN-IQ PN-DOMP sweep"; "Sparse PNNN N12"];
rows = sweep.results([],:);
for index = 1:numel(expectedModels)
    match = sweep.results.TargetRealParameters == target & ...
        string(sweep.results.Model) == expectedModels(index);
    if nnz(match) ~= 1
        error('run_fair_PNNN_vs_PNGMP_DOMP:MissingSweepPoint', ...
            'The selected target must contain exactly the three sweep families.');
    end
    rows = [rows; sweep.results(match,:)]; %#ok<AGROW>
end
allTargetRows = sweep.results.TargetRealParameters == target;
if nnz(allTargetRows) ~= 3 || ...
        any(string(rows.SweepRole) ~= "Sweep point") || ...
        any(rows.ActualRealParameters ~= target)
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidSweepPoint', ...
        'The selected target must contain three exact-parameter sweep rows.');
end
pnnnRow = string(rows.Model) == "Sparse PNNN N12";
if rows.ActiveWeights(pnnnRow) + rows.ActiveBiases(pnnnRow) ~= target
    error('run_fair_PNNN_vs_PNGMP_DOMP:PNNNParameterMismatch', ...
        'Sparse PNNN weights and biases must equal the selected target.');
end
end

function artifact = loadSignedArtifact(filename, identity, signature)
filename = normalizePath(filename);
if ~isfile(filename)
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingArtifact', ...
        'Required sweep artifact not found: %s', filename);
end
variables = whos('-file', filename);
if ~any(strcmp({variables.name}, 'checkpointArtifact'))
    error('run_fair_PNNN_vs_PNGMP_DOMP:UnsignedArtifact', ...
        'Sweep artifact does not contain checkpointArtifact: %s', filename);
end
loaded = load(filename, 'checkpointArtifact');
artifact = loaded.checkpointArtifact;
required = {'schemaVersion','sweepIdentity','experimentSignature','payload'};
if ~isstruct(artifact) || ~all(isfield(artifact, required)) || ...
        artifact.schemaVersion ~= 1 || ...
        ~isequaln(artifact.sweepIdentity, identity) || ...
        ~isequaln(artifact.experimentSignature, signature)
    error('run_fair_PNNN_vs_PNGMP_DOMP:IncompatibleArtifact', ...
        'Sweep artifact identity or experiment signature is incompatible.');
end
end

function validateSummary(payload, selectedRows, target)
if ~isstruct(payload) || ~isfield(payload, 'results') || ...
        ~istable(payload.results)
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidSummary', ...
        'The sweep summary does not contain its canonical result table.');
end
summaryRows = payload.results(payload.results.TargetRealParameters == target,:);
columns = {'Model','ActualRealParameters','FullSignalNMSEdB','FLOPsPerSample'};
if height(summaryRows) ~= 3 || ...
        ~isequaln(summaryRows(:,columns), selectedRows(:,columns))
    error('run_fair_PNNN_vs_PNGMP_DOMP:SummaryMismatch', ...
        'The supplied sweep and signed summary disagree at the selected target.');
end
end

function [predictions, supports] = selectLinearPoint(payload, rows, target)
required = {'complexTable','pnTable','supports','paths','predictions'};
if ~isstruct(payload) || ~all(isfield(payload, required)) || ...
        ~all(isfield(payload.predictions, {'complexFull','pnFull'}))
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidLinearArtifact', ...
        'The linear artifact lacks reusable sweep predictions.');
end
complexIndex = uniqueIndex(payload.complexTable, target);
pnIndex = uniqueIndex(payload.pnTable, target);
assertRowMatches(payload.complexTable(complexIndex,:), rows(1,:));
assertRowMatches(payload.pnTable(pnIndex,:), rows(2,:));
predictions = struct( ...
    'complexGMP', payload.predictions.complexFull(:,complexIndex), ...
    'pnIQ', payload.predictions.pnFull(:,pnIndex));
supports = struct( ...
    'complex', payload.supports.complex{complexIndex}, ...
    'pnFeatures', payload.supports.pnFeatures{pnIndex}, ...
    'pnComplex', payload.supports.pnComplex{pnIndex});
end

function prediction = selectPNNNPoint(payload, densePayload, rows, target)
required = {'target','row','fullSignalPrediction','denseSourceSignature'};
if ~isstruct(payload) || ~all(isfield(payload, required)) || ...
        payload.target ~= target || ~istable(payload.row) || ...
        height(payload.row) ~= 1 || ~isstruct(densePayload) || ...
        ~isfield(densePayload, 'signature') || ...
        ~isequaln(payload.denseSourceSignature, densePayload.signature)
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidPNNNArtifact', ...
        'The sparse PNNN artifact is not compatible with its dense source.');
end
assertRowMatches(payload.row, rows(3,:));
prediction = payload.fullSignalPrediction(:);
end

function index = uniqueIndex(tableValue, target)
match = tableValue.TargetRealParameters == target;
if nnz(match) ~= 1
    error('run_fair_PNNN_vs_PNGMP_DOMP:MissingArtifactPoint', ...
        'The selected target is not unique in a sweep artifact.');
end
index = find(match);
end

function assertRowMatches(actual, expected)
columns = {'Model','ActualRealParameters','FullSignalNMSEdB','FLOPsPerSample'};
if ~isequaln(actual(:,columns), expected(:,columns))
    error('run_fair_PNNN_vs_PNGMP_DOMP:ArtifactRowMismatch', ...
        'A checkpoint row disagrees with the signed sweep summary.');
end
end

function validatePredictionSizes( ...
    linearPredictions, pnnnPrediction, target, indices)
count = numel(linearPredictions.complexGMP);
valid = count > 0 && numel(linearPredictions.pnIQ) == count && ...
    numel(pnnnPrediction) == count && numel(target) == count && ...
    numel(indices) == count && ...
    all(isfinite(linearPredictions.complexGMP)) && ...
    all(isfinite(linearPredictions.pnIQ)) && all(isfinite(pnnnPrediction));
if ~valid
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidPredictions', ...
        'The three stored full-signal predictions must be aligned and finite.');
end
end

function validateFixedPredictionSizes(tableValue, predictions, count)
models = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"];
lambdas = [1e-3; 1e-4; 1e-5];
families = {'complexGMP','pnIQ'};
fields = {'lambda1e3','lambda1e4','lambda1e5'};
valid = istable(tableValue) && height(tableValue) == 6 && ...
    all(ismember({'Model','FixedLambda'}, ...
    tableValue.Properties.VariableNames)) && ...
    isequal(sort(unique(string(tableValue.Model))), sort(models)) && ...
    isequal(sort(unique(tableValue.FixedLambda)), sort(lambdas));
for familyIndex = 1:numel(families)
    family = families{familyIndex};
    valid = valid && isfield(predictions, family) && ...
        all(isfield(predictions.(family), fields));
    if valid
        for fieldIndex = 1:numel(fields)
            value = predictions.(family).(fields{fieldIndex});
            valid = valid && numel(value) == count && all(isfinite(value));
        end
    end
end
if ~valid
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidFixedPredictions', ...
        'Exactly six aligned fixed-lambda predictions are required.');
end
end

function value = normalizePath(value)
valid = (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value));
if ~valid || ismissing(string(value)) || strlength(string(value)) == 0
    error('run_fair_PNNN_vs_PNGMP_DOMP:InvalidPath', ...
        'Sweep artifact paths must be nonempty text scalars.');
end
value = char(string(value));
end

function valid = validSignature(value)
valid = isstruct(value) && isfield(value, 'schemaVersion') && ...
    isfield(value, 'digest') && isscalar(string(value.digest)) && ...
    strlength(string(value.digest)) > 0;
end
