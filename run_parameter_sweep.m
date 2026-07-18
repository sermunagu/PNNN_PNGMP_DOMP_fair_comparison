function result = run_parameter_sweep(parameterGrid, resultsRoot)
% run_parameter_sweep - Run or resume the three-family complexity sweep.
% The scientific stages appear in execution order, while signed checkpoints
% keep linear, fixed-ridge, dense N12, and sparse targets independent.

%% Configure the experiment
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'config'));
for folder = ["complexity","domp","metrics","pn_gmp_comparison", ...
        "pnnn","pnnn/pruning","splits","sweep"]
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

if ~isfile(linearFile)
    fprintf('[Linear] Building designs, DOMP paths, fits, and predictions...\n');
    
    linear = run_linear_sweep(x, y, split, cfg);
    checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
        'payload', linear);
    
    save(linearFile, 'checkpointArtifact', '-v7.3');
    fprintf('[Linear] Completed.\n');

else
    saved = load(linearFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;
    
    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', ...
            'The linear checkpoint belongs to another experiment.');
    end
    linear = checkpointArtifact.payload;
    
    if any(string(linear.complexTable.Model) ~= ...
            "Complex GMP DOMP sweep") || ...
            any(string(linear.pnTable.Model) ~= ...
            "Independent PN-IQ PN-DOMP sweep")
        error('run_parameter_sweep:FamilyMismatch', 'The linear checkpoint mixes model families.');
    end
    fprintf('[Linear] Reused matrices, paths, fits, and predictions.\n');
end

%% Fixed Ridge: refit both linear families on the stored supports
fixedFile = fullfile(resultDirectory, 'fixed_lambda_linear_sweep.mat');

if ~isfile(fixedFile)
    fixedLinear = run_fixed_ridge_sweep(x, y, split, cfg, linear);
    fixedResults = fixedLinear.table;
    checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
        'payload', fixedResults);
    save(fixedFile, 'checkpointArtifact', '-v7.3');
    fprintf('[Fixed ridge] Completed.\n');

else
    saved = load(fixedFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;
    
    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', ...
            'The fixed-ridge checkpoint belongs to another experiment.');
    end
    fixedResults = checkpointArtifact.payload;
    
    if ~isequal(sort(unique(string(fixedResults.Model))), ...
            sort(["Complex GMP-DOMP"; "PN-IQ PN-DOMP"]))
        error('run_parameter_sweep:FamilyMismatch', ...
            'The fixed-ridge checkpoint mixes model families.');
    end
    
    fprintf('[Fixed ridge] Reused supplementary checkpoint.\n');
end

%% Sparse PNNN: prepare or reuse one immutable dense N12 source
denseFile = fullfile(resultDirectory, 'sweep_dense_source.mat');

if ~isfile(denseFile)
    fprintf('[PNNN] Selecting epochs and fitting one dense N12 source...\n');
    [denseSource, features, neuralTargets, phaseRotation] = ...
        prepare_pnnn_dense_source( ...
            x, y, split, cfg, cfg.reducedRealParameterTarget);

    checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
        'payload', denseSource);
    save(denseFile, 'checkpointArtifact', '-v7.3');
    fprintf('[PNNN] Dense source completed.\n');

else
    saved = load(denseFile, 'checkpointArtifact');
    checkpointArtifact = saved.checkpointArtifact;
    
    if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
        error('run_parameter_sweep:ExperimentMismatch', ...
            'The dense PNNN checkpoint belongs to another experiment.');
    end
    
    denseSource = checkpointArtifact.payload;
    [features, neuralTargets, phaseRotation] = buildPhaseNormDataset( ...
        x, y, cfg.pnnn.M, cfg.pnnn.orders, cfg.pnnn.featMode);
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
        point = fit_sparse_pnnn_target(denseSource, target, features, ...
            neuralTargets, phaseRotation, y, split, cfg);
        checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
            'payload', point);
        save(filename, 'checkpointArtifact', '-v7.3');
        fprintf('[PNNN %d/%d] Completed.\n', index, numel(targets));
    
    else
        saved = load(filename, 'checkpointArtifact');
        checkpointArtifact = saved.checkpointArtifact;
        
        if ~isequaln(checkpointArtifact.sweepIdentity, sweepIdentity)
            error('run_parameter_sweep:ExperimentMismatch', ...
                'The sparse PNNN checkpoint belongs to another experiment.');
        end
        
        if ~ismember(target, checkpointArtifact.sweepIdentity.parameterGrid)
            error('run_parameter_sweep:UnsignedTarget', ...
                'The requested target is not part of the signed grid.');
        end
        
        point = checkpointArtifact.payload;
        
        if point.row.TargetRealParameters ~= target || ...
                string(point.row.Model) ~= "Sparse PNNN N12"
            error('run_parameter_sweep:FamilyMismatch', ...
                'The sparse checkpoint belongs to another target or family.');
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

%% Write the canonical CSV files and two sweep figures
fprintf('[Output] Writing complexity_sweep.csv...\n');
writeTableAtomically(results, resultDirectory, 'complexity_sweep.csv');
fprintf('[Output] Writing fixed_lambda_linear_sweep.csv...\n');
writeTableAtomically(fixedResults, resultDirectory, ...
    'fixed_lambda_linear_sweep.csv');
fprintf('[Output] Writing NMSE-vs-parameters figure...\n');
plotSweepFigure(results, fixedResults, 'ActualRealParameters', ...
    'Active real parameters', fullfile(resultDirectory, ...
    'comparison_nmse_parameters_sweep.png'), true, ...
    cfg.fixedRidgeLambdas);
fprintf('[Output] Writing NMSE-vs-FLOPs figure...\n');
plotSweepFigure(results, table(), 'FLOPsPerSample', ...
    'FLOPs per sample', fullfile(resultDirectory, ...
    'comparison_nmse_flops_sweep.png'), false, ...
    cfg.fixedRidgeLambdas);

%% Save the signed scientific table without duplicating checkpoint contents
summaryFile = fullfile(resultDirectory, 'complexity_sweep.mat');
summaryPayload = struct('results', results, ...
    'evaluationProtocol', evaluationProtocol);
checkpointArtifact = struct('sweepIdentity', sweepIdentity, ...
    'payload', summaryPayload);
save(summaryFile, 'checkpointArtifact', '-v7.3');
result = struct('results', results, 'sweepIdentity', sweepIdentity, ...
    'evaluationProtocol', evaluationProtocol, ...
    'resultDirectory', string(resultDirectory));
disp(results(:, {'Model','ActualRealParameters', ...
    'SelectedLambda','FullSignalNMSEdB','FLOPsPerSample'}));
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

function plotSweepFigure(results, fixed, xVariable, xLabel, filename, ...
    includeFixed, fixedLambdas)
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 900 600]);
cleanup = onCleanup(@() close(figureHandle));
hold on;
models = ["Complex GMP DOMP sweep", ...
    "Independent PN-IQ PN-DOMP sweep", "Sparse PNNN N12"];
labels = ["Complex GMP-DOMP", "PN-IQ PN-DOMP", "Sparse PNNN N12"];
mainLines = gobjects(3, 1);

for index = 1:3
    rows = string(results.Model) == models(index);
    x = results.(xVariable)(rows);
    y = results.FullSignalNMSEdB(rows);
    [x, order] = sort(x);
    mainLines(index) = plot(x, y(order), '-o', ...
        'LineWidth', 1 + includeFixed, 'DisplayName', labels(index));
end

if includeFixed
    fixedModels = ["Complex GMP-DOMP", "PN-IQ PN-DOMP"];
    styles = ["--", ":", "-."];
    for modelIndex = 1:2
        for lambdaIndex = 1:numel(fixedLambdas)
            rows = string(fixed.Model) == fixedModels(modelIndex) & ...
                fixed.FixedLambda == fixedLambdas(lambdaIndex);
            x = fixed.(xVariable)(rows);
            y = fixed.FullSignalNMSEdB(rows);
            [x, order] = sort(x);
            plot(x, y(order), 'LineStyle', styles(lambdaIndex), ...
                'Color', mainLines(modelIndex).Color, 'LineWidth', 1, ...
                'DisplayName', fixedModels(modelIndex) + ", lambda=" + ...
                compose('%g', fixedLambdas(lambdaIndex)));
        end
    end
end

grid on;
xlabel(xLabel);
ylabel('Full-signal NMSE (dB)');
legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
    'NumColumns', 3, 'Interpreter', 'none');
exportgraphics(figureHandle, filename, 'Resolution', 160, ...
    'BackgroundColor', 'white');
clear cleanup;
end

function deleteIfPresent(filename)
if isfile(filename)
    delete(filename);
end
end
