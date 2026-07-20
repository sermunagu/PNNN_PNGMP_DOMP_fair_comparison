function result = run_parameter_sweep(parameterGrid, resultsRoot)
% run_parameter_sweep - Run or resume the three-family complexity sweep.
% The scientific stages appear in execution order, while signed checkpoints
% keep linear, fixed-ridge, dense N12, and sparse targets independent.

%% Configure the experiment
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
for folder = ["complexity","domp","metrics","pn_gmp_comparison", "pnnn","pnnn/pruning","splits","sweep"]
    addpath(fullfile(projectRoot, 'toolbox', folder));
end

cfg = getFairDOMPComparisonConfig(projectRoot);
if nargin >= 2 && ~isempty(resultsRoot)
    cfg.sweep.resultsRoot = char(string(resultsRoot));
end

if nargin >= 1 && ~isempty(parameterGrid)
    cfg.sweep.parameterGrid = unique(double(parameterGrid(:).'));
end

targets = cfg.sweep.parameterGrid;

%% Verify the paper-export toolchain before data loading or training
fprintf('[Preflight] Checking matlab2tikz and standalone LaTeX export...\n');
paperToolchain = preflightPaperFigureToolchain(projectRoot, cfg.paper);
fprintf('[Preflight] Paper-figure toolchain passed.\n');

%% Load the measurement and build the common split
measurement = load(cfg.measurementFile, 'x', 'y');
[x, y] = selectXYByMapping(measurement.x, measurement.y, cfg.mappingMode);
x = x(:);
y = y(:);

experimentSignature = buildExperimentSignature(x, y, cfg);

if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end

split = buildCommonComparisonSplit(x, y, cfg);

sweepIdentity = buildSweepIdentity(cfg, experimentSignature);
resultDirectory = resolveResultDirectory(cfg, sweepIdentity);

if ~isfolder(resultDirectory)
    mkdir(resultDirectory);
end

fprintf('\n=== Signed parameter-complexity sweep ===\n');
fprintf('Result directory: %s\n', resultDirectory);

%% Complex GMP and PN-IQ: select paths, fit prefixes, and predict
linearFile = fullfile(resultDirectory, 'linear_sweep.mat');
linearNeedsRebuild = ~isfile(linearFile);
linearWasRebuilt = false;

if ~linearNeedsRebuild
    saved = load(linearFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;

    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', 'The linear checkpoint belongs to another experiment.');
    end
    linear = checkpointArtifact.payload;

    if ~hasLinearCoefficientDefinition(linear, ...
            cfg.sweep.coefficientRangeDefinition)
        linearNeedsRebuild = true;
        fprintf(['[Linear] Checkpoint coefficient-range definition is ' ...
            'missing or stale; rebuilding linear models only.\n']);
    else
        if any(string(linear.complexTable.Model) ~= "Complex GMP DOMP sweep") || ...
                any(string(linear.pnTable.Model) ~= "Independent PN-IQ PN-DOMP sweep")
            error('run_parameter_sweep:FamilyMismatch', 'The linear checkpoint mixes model families.');
        end
        fprintf('[Linear] Reused matrices, paths, fits, and predictions.\n');
    end
end

if linearNeedsRebuild
    fprintf('[Linear] Building designs, DOMP paths, fits, and predictions...\n');
    linear = run_linear_sweep(x, y, split, cfg);
    checkpointArtifact = struct('sweepIdentity', sweepIdentity, 'payload', linear);
    save(linearFile, 'checkpointArtifact', '-v7.3');
    linearWasRebuilt = true;
    fprintf('[Linear] Completed.\n');
end

%% Fixed Ridge: refit both linear families on the stored supports
fixedFile = fullfile(resultDirectory, 'fixed_lambda_linear_sweep.mat');
fixedNeedsRebuild = linearWasRebuilt || ~isfile(fixedFile);

if ~fixedNeedsRebuild
    saved = load(fixedFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;

    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', 'The fixed-ridge checkpoint belongs to another experiment.');
    end
    fixedResults = checkpointArtifact.payload;

    if ~hasFixedCoefficientDefinition(fixedResults, ...
            cfg.sweep.coefficientRangeDefinition)
        fixedNeedsRebuild = true;
        fprintf(['[Fixed ridge] Checkpoint coefficient-range definition is ' ...
            'missing or stale; rebuilding fixed Ridge only.\n']);
    else
        if ~isequal(sort(unique(string(fixedResults.Model))), sort(["Complex GMP-DOMP"; "PN-IQ PN-DOMP"]))
            error('run_parameter_sweep:FamilyMismatch', 'The fixed-ridge checkpoint mixes model families.');
        end
        fprintf('[Fixed ridge] Reused supplementary checkpoint.\n');
    end
end

if fixedNeedsRebuild
    fixedLinear = run_fixed_ridge_sweep(x, y, split, cfg, linear);
    fixedResults = fixedLinear.table;
    checkpointArtifact = struct('sweepIdentity', sweepIdentity, 'payload', fixedResults);
    save(fixedFile, 'checkpointArtifact', '-v7.3');
    fprintf('[Fixed ridge] Completed.\n');
end

%% Sparse PNNN: prepare or reuse one immutable dense N12 source
denseFile = fullfile(resultDirectory, 'sweep_dense_source.mat');

if ~isfile(denseFile)
    fprintf('[PNNN] Selecting epochs and fitting one dense N12 source...\n');
    [denseSource, features, neuralTargets, phaseRotation] = prepare_pnnn_dense_source( ...
            x, y, split, cfg, cfg.reducedRealParameterTarget);

    checkpointArtifact = struct('sweepIdentity', sweepIdentity, 'payload', denseSource);
    save(denseFile, 'checkpointArtifact', '-v7.3');
    fprintf('[PNNN] Dense source completed.\n');

else
    saved = load(denseFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;
    
    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', 'The dense PNNN checkpoint belongs to another experiment.');
    end
    
    denseSource = checkpointArtifact.payload;
    [features, neuralTargets, phaseRotation] = buildPhaseNormDataset(x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
    features = features.';
    neuralTargets = neuralTargets.';
    phaseRotation = phaseRotation(:);
    fprintf('[PNNN] Reused dense N12 checkpoint.\n');
end

% Console verbosity belongs to this execution, not to the dense checkpoint.
denseSource.runtimeConfig.training.verbose = cfg.training.verbose;

%% Sparse PNNN: prune and predict every requested parameter budget
pnnnRows = linear.complexTable([],:);

for index = 1:numel(targets)
    target = targets(index);
    fprintf('[PNNN %d/%d] Target %d parameters...\n', index, numel(targets), target);
    
    artifactName = compose("pnnn_target_%04d.mat", target);
    filename = fullfile(resultDirectory, artifactName);
    
    if ~isfile(filename)
        point = fit_sparse_pnnn_target(denseSource, target, features, neuralTargets, phaseRotation, y, split, cfg);
        checkpointArtifact = struct('sweepIdentity', sweepIdentity, 'payload', point);
        save(filename, 'checkpointArtifact', '-v7.3');
        fprintf('[PNNN %d/%d] Completed.\n', index, numel(targets));
    
    else
        saved = load(filename, 'checkpointArtifact');
        checkpointArtifact = saved.checkpointArtifact;
        
        if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
            error('run_parameter_sweep:ExperimentMismatch', 'The sparse PNNN checkpoint belongs to another experiment.');
        end
        
        if ~ismember(target, checkpointArtifact.sweepIdentity.parameterGrid)
            error('run_parameter_sweep:UnsignedTarget', 'The requested target is not part of the signed grid.');
        end
        
        point = checkpointArtifact.payload;
        
        if point.row.TargetRealParameters ~= target || string(point.row.Model) ~= "Sparse PNNN N12"
            error('run_parameter_sweep:FamilyMismatch', 'The sparse checkpoint belongs to another target or family.');
        end
        
        fprintf('[PNNN %d/%d] Reused checkpoint.\n', index, numel(targets));
    end
    
    pnnnRows(index, :) = point.row;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES AND METADATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Build the two public tables
results = [linear.complexTable; linear.pnTable; pnnnRows];
evaluationProtocol = struct( ...
    'internalValidationUsedForSelectionOnly', true, ...
    'completeSignalIncludesIdentification', ...
        all(ismember(split.identificationIndices, split.fullSignalIndices)), ...
    'completeSignalIsIndependentHoldout', ...
        ~any(ismember(split.identificationIndices, split.fullSignalIndices)));
expectedModels = sort(["Complex GMP DOMP sweep"; ...
    "Independent PN-IQ PN-DOMP sweep"; "Sparse PNNN N12"]);
if height(results) ~= 3*numel(targets) || ...
        height(fixedResults) ~= 6*numel(targets)
    error('run_parameter_sweep:ResultRowCount', ...
        'The principal and fixed-ridge tables must contain 3N and 6N rows.');
end
for target = targets
    rows = results.TargetRealParameters == target;
    if nnz(rows) ~= 3 || ...
            ~isequal(sort(string(results.Model(rows))), expectedModels) || ...
            any(results.ActualRealParameters(rows) ~= target)
        error('run_parameter_sweep:ResultTargetMismatch', ...
            'Each target must contain the three exact-parameter families.');
    end
end
pnnnMask = string(results.Model) == "Sparse PNNN N12";
if any(results.ActiveWeights(pnnnMask) + results.ActiveBiases(pnnnMask) ~= ...
        results.ActualRealParameters(pnnnMask))
    error('run_parameter_sweep:PNNNParameterMismatch', ...
        'PNNN active weights and biases must equal its parameter count.');
end
if width(results) ~= 13 || width(fixedResults) ~= 8 || ...
        any(~isfinite(results.MaxAbsRealParameter)) || ...
        any(~isfinite(fixedResults.MaxAbsRealParameter))
    error('run_parameter_sweep:CoefficientRangeContract', ...
        ['The schema-v3 principal/fixed tables must contain 13/8 columns ' ...
        'and finite MaxAbsRealParameter values.']);
end

%% Select the first jointly stabilized minimum-complexity point
selection = selectOperatingPoint(results, cfg.selection);
coefficientMetadata = struct( ...
    'definition', string(cfg.sweep.coefficientRangeDefinition), ...
    'description', "Maximum absolute active real scalar in normalized " + ...
        "input/output parameterization. Linear-family coefficients are " + ...
        "expressed for unit-peak identification input/output and unit-2-norm " + ...
        "regressor columns; " + ...
        "PNNN parameters remain in their z-score-normalized parameterization.", ...
    'warning', "The metric compares numerical dynamic range for " + ...
        "implementation and quantization. Linear families use unit-peak " + ...
        "identification I/O and unit-norm regressors, whereas PNNN uses its " + ...
        "z-score parameterization; magnitudes must not be interpreted as " + ...
        "identical physical gains.");
warning('run_parameter_sweep:ParameterizationWarning', '%s', ...
    coefficientMetadata.warning);
fprintf('%s\n', selection.summarySentence);

%% Write canonical tables, selection diagnostics, and paper figures
fprintf('[Output] Writing complexity_sweep.csv...\n');
writeTableAtomically(results, resultDirectory, 'complexity_sweep.csv');
fprintf('[Output] Writing fixed_lambda_linear_sweep.csv...\n');
writeTableAtomically(fixedResults, resultDirectory, ...
    'fixed_lambda_linear_sweep.csv');
fprintf('[Output] Writing operating-point selection diagnostics...\n');
writeTableAtomically(selection.diagnosticsTable, resultDirectory, ...
    'operating_point_selection.csv');
writeTableAtomically(selection.sensitivityTable, resultDirectory, ...
    'operating_point_selection_sensitivity.csv');
writeTextAtomically([selection.summarySentence; ...
    coefficientMetadata.warning], resultDirectory, ...
    'operating_point_selection_summary.txt');
exportOptions = struct('latexmkCommand', cfg.paper.latexmkCommand);
nmseOptions = struct('metricVariable', 'FullSignalNMSEdB', ...
    'metricLabel', 'Full-signal NMSE (dB)', 'includeFixed', true, ...
    'fixedLambdas', cfg.fixedRidgeLambdas, 'isNMSE', true, ...
    'annotateSelected', false, 'selection', selection, ...
    'exportOptions', exportOptions);
fprintf('[Output] Writing NMSE-vs-parameters figure...\n');
figureFiles.nmseParameters = plotSweepPaperFigure(results, fixedResults, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(resultDirectory, 'comparison_nmse_parameters_sweep'), ...
    nmseOptions);
fprintf('[Output] Writing NMSE-vs-FLOPs figure...\n');
nmseOptions.includeFixed = false;
figureFiles.nmseFLOPs = plotSweepPaperFigure(results, table(), ...
    'FLOPsPerSample', 'FLOPs per sample', ...
    fullfile(resultDirectory, 'comparison_nmse_flops_sweep'), ...
    nmseOptions);
fprintf('[Output] Writing coefficient-range figure...\n');
rangeOptions = struct('metricVariable', 'MaxAbsRealParameter', ...
    'metricLabel', ['Max. abs. real coefficient ' ...
        '(unit-peak I/O; unit-norm regressors)'], ...
    'includeFixed', true, 'fixedLambdas', cfg.fixedRidgeLambdas, ...
    'useLogWhenPositive', true, 'annotateSelected', false, ...
    'exportOptions', exportOptions);
figureFiles.maxAbsParameter = plotSweepPaperFigure(results, fixedResults, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(resultDirectory, 'comparison_max_abs_parameter_sweep'), ...
    rangeOptions);

%% Save the signed scientific table without duplicating checkpoint contents
summaryFile = fullfile(resultDirectory, 'complexity_sweep.mat');
summaryPayload = struct('results', results, ...
    'evaluationProtocol', evaluationProtocol, ...
    'selection', selection, ...
    'coefficientMetadata', coefficientMetadata, ...
    'figureFiles', figureFiles);
checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
    'payload', summaryPayload);
save(summaryFile, 'checkpointArtifact', '-v7.3');
result = struct('results', results, 'sweepIdentity', sweepIdentity, ...
    'evaluationProtocol', evaluationProtocol, ...
    'selection', selection, ...
    'coefficientMetadata', coefficientMetadata, ...
    'figureFiles', figureFiles, ...
    'paperToolchain', paperToolchain, ...
    'selectedParameters', selection.selectedParameters, ...
    'criterionName', selection.criterionName, ...
    'stabilizationWindowParameters', ...
        selection.stabilizationWindowParameters, ...
    'stabilizationToleranceDb', selection.stabilizationToleranceDb, ...
    'selectedComplexGMPFutureGainDb', ...
        selection.selectedComplexGMPFutureGainDb, ...
    'selectedPNIQFutureGainDb', selection.selectedPNIQFutureGainDb, ...
    'selectedSparsePNNNFutureGainDb', ...
        selection.selectedSparsePNNNFutureGainDb, ...
    'selectedWorstFutureGainDb', selection.selectedWorstFutureGainDb, ...
    'selectedComplexGMPFLOPs', selection.selectedComplexGMPFLOPs, ...
    'selectedPNIQFLOPs', selection.selectedPNIQFLOPs, ...
    'selectedSparsePNNNFLOPs', selection.selectedSparsePNNNFLOPs, ...
    'resultDirectory', string(resultDirectory));
disp(results(:, {'Model','ActualRealParameters', ...
    'SelectedLambda','FullSignalNMSEdB','FLOPsPerSample', ...
    'MaxAbsRealParameter'}));
fprintf('Sweep completed.\n');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUXILIAR FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function identity = buildSweepIdentity(cfg, experimentSignature)
trainingIdentity = struct( ...
    'optimizer', cfg.training.optimizer, ...
    'miniBatchSize', cfg.training.miniBatchSize, ...
    'initialLearnRate', cfg.training.initialLearnRate, ...
    'learnRateSchedule', cfg.training.learnRateSchedule, ...
    'learnRateDropFactor', cfg.training.learnRateDropFactor, ...
    'shuffle', cfg.training.shuffle, ...
    'executionEnvironment', cfg.training.executionEnvironment, ...
    'inputDataFormats', cfg.training.inputDataFormats, ...
    'targetDataFormats', cfg.training.targetDataFormats, ...
    'historicalTrainFraction', cfg.training.historicalTrainFraction, ...
    'historicalMaxEpochs', cfg.training.historicalMaxEpochs, ...
    'historicalLearnRateDropPeriod', ...
        cfg.training.historicalLearnRateDropPeriod, ...
    'historicalValidationPatience', ...
        cfg.training.historicalValidationPatience);
pruningIdentity = struct( ...
    'historicalFineTuneEpochs', cfg.pruning.historicalFineTuneEpochs, ...
    'fineTuneInitialLearnRate', cfg.pruning.fineTuneInitialLearnRate, ...
    'fineTuneSeedOffset', cfg.pruning.fineTuneSeedOffset);
identity = struct('schemaVersion', cfg.sweep.schemaVersion, ...
    'experimentSignature', experimentSignature, ...
    'parameterGrid', cfg.sweep.parameterGrid, ...
    'lambdaGrid', cfg.lambdaGrid, ...
    'fixedRidgeLambdas', cfg.fixedRidgeLambdas, ...
    'gmpPopulation', ...
        [cfg.gmp.Qpmax cfg.gmp.Qnmax cfg.gmp.Pmax cfg.gmp.maxPopulation], ...
    'dompOptions', cfg.gmp.dompOptions, ...
    'denseSelectionTarget', cfg.reducedRealParameterTarget, ...
    'pnnn', cfg.pnnn, 'training', trainingIdentity, ...
    'pruning', pruningIdentity);

engine = javaMethod('getInstance', ...
    'java.security.MessageDigest', 'SHA-256');

bytes = unicode2native(jsonencode(identity), 'UTF-8');
engine.update(typecast(uint8(bytes), 'int8'));
digest = typecast(int8(engine.digest()), 'uint8');
identity.digest = string(lower(reshape(dec2hex(digest, 2).', 1, [])));
end

function directory = resolveResultDirectory(cfg, identity)
name = "sweep_" + extractBefore(string(identity.digest), 13);

if ~cfg.sweep.resume
    name = name + "_" + string(datetime('now', ...
        'Format', 'yyyyMMdd_HHmmss'));
end
directory = fullfile(cfg.sweep.resultsRoot, name);
end

function writeTableAtomically(value, directory, filename)
finalFile = fullfile(directory, filename);
temporaryFile = [tempname(directory) '.csv'];
cleanup = onCleanup(@() deleteIfPresent(temporaryFile));
writetable(value, temporaryFile);
[moved, message] = movefile(temporaryFile, finalFile, 'f');

if ~moved
    error('run_parameter_sweep:CSVMoveFailed', ...
        'Could not install %s: %s', filename, message);
end

clear cleanup;
end

function writeTextAtomically(value, directory, filename)
finalFile = fullfile(directory, filename);
temporaryFile = [tempname(directory) '.txt'];
cleanup = onCleanup(@() deleteIfPresent(temporaryFile));
writelines(string(value), temporaryFile);
[moved, message] = movefile(temporaryFile, finalFile, 'f');
if ~moved
    error('run_parameter_sweep:TextMoveFailed', ...
        'Could not install %s: %s', filename, message);
end
clear cleanup;
end

function deleteIfPresent(filename)
if isfile(filename)
    delete(filename);
end
end

function matches = hasLinearCoefficientDefinition(payload, expected)
matches = isstruct(payload) && ...
    isfield(payload, 'coefficientRangeDefinition');
if matches
    definition = string(payload.coefficientRangeDefinition);
    matches = isscalar(definition) && definition == string(expected);
end
end

function matches = hasFixedCoefficientDefinition(value, expected)
matches = istable(value) && isstruct(value.Properties.UserData) && ...
    isfield(value.Properties.UserData, 'coefficientRangeDefinition');
if matches
    definition = string( ...
        value.Properties.UserData.coefficientRangeDefinition);
    matches = isscalar(definition) && definition == string(expected);
end
end
