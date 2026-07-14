function [pruningState, stats] = createMagnitudePruningMasks(net, pruningCfg)
% createMagnitudePruningMasks - Build magnitude pruning masks.
%
% This function ranks prunable learnable parameters by absolute value and
% creates binary masks for PNNN fine-tuning. It supports global and layerwise
% scopes, and targets by sparsity or final active trainable-parameter count.
%
% Inputs:
%   net - Trained dlnetwork returned by trainnet.
%   pruningCfg - Validated pruning configuration struct.
%
% Outputs:
%   pruningState - Struct containing per-learnable binary masks.
%   stats - Pruning statistics for metadata.

if ~isa(net, 'dlnetwork')
    error("Magnitude pruning con mascaras requiere que trainnet devuelva un dlnetwork.");
end

pruningCfg = normalizePruningCfgForMasks(pruningCfg);
if string(pruningCfg.structureMode) ~= "unstructured"
    [pruningState, stats] = createStructuredPruningMasks(net, pruningCfg);
    return;
end
learnables = net.Learnables;
nLearnables = height(learnables);
masks = cell(nLearnables, 1);
candidates = struct("row", {}, "numel", {}, "size", {}, "name", {}, ...
    "parameter", {}, "magnitudes", {});

totalTrainableParams = 0;
protectedTrainableParams = 0;

for i = 1:nLearnables
    value = learnables.Value{i};
    data = learnableToNumeric(value);
    masks{i} = true(size(data));
    paramCount = numel(data);
    totalTrainableParams = totalTrainableParams + paramCount;

    parameterName = string(learnables.Parameter(i));
    if isPrunableParameter(parameterName, pruningCfg.includeBiases)
        candidates(end+1).row = i; %#ok<AGROW>
        candidates(end).numel = paramCount;
        candidates(end).size = size(data);
        candidates(end).name = learnableName(learnables, i);
        candidates(end).parameter = lower(parameterName);
        candidates(end).magnitudes = abs(data(:));
    else
        protectedTrainableParams = protectedTrainableParams + paramCount;
    end
end

stats = initPruningStats(pruningCfg);
totalPrunableParams = sum([candidates.numel]);
if totalPrunableParams == 0
    error("No se encontraron parametros podables para pruning.");
end

[targetActiveTrainableParams, targetActivePrunableParams, ...
    prunedPrunableParams, targetPrunableSparsity] = resolveTargetCounts( ...
    pruningCfg, totalTrainableParams, totalPrunableParams, ...
    protectedTrainableParams);

scope = string(pruningCfg.scope);
if scope == "layerwise"
    [masks, parameterNames, parameterTotal, parameterPruned] = ...
        createLayerwiseMasks(masks, candidates, prunedPrunableParams);
else
    [masks, parameterNames, parameterTotal, parameterPruned] = ...
        createGlobalMasks(masks, candidates, prunedPrunableParams);
end

actualPrunedPrunableParams = sum(parameterPruned);
remainingPrunableParams = totalPrunableParams - actualPrunedPrunableParams;
paramCounts = summarizeTrainableParameters(net, pruningCfg.includeBiases, masks);
activeWeightParams = paramCounts.activeWeightParams;
activeBiasParams = paramCounts.activeBiasParams;
actualActiveTrainableParams = paramCounts.actualActiveTrainableParams;
actualPrunableSparsity = paramCounts.actualPrunableSparsity;

stats.scope = char(scope);
stats.includeBiases = pruningCfg.includeBiases;
stats.targetMode = char(pruningCfg.targetMode);
stats.structureMode = char(pruningCfg.structureMode);
stats.structuredRanking = char(pruningCfg.structuredRanking);
stats.structuredTargetPolicy = char(pruningCfg.structuredTargetPolicy);
stats.hybridExactTarget = logical(pruningCfg.hybridExactTarget);
stats.sparsityTarget = targetPrunableSparsity;
stats.sparsityActual = actualPrunableSparsity;
stats.totalPodableParams = totalPrunableParams;
stats.totalTrainableParams = totalTrainableParams;
stats.totalPrunableParams = totalPrunableParams;
stats.protectedTrainableParams = protectedTrainableParams;
stats.totalWeightParams = paramCounts.totalWeightParams;
stats.totalBiasParams = paramCounts.totalBiasParams;
stats.targetActiveTrainableParams = targetActiveTrainableParams;
stats.targetActivePrunableParams = targetActivePrunableParams;
stats.targetActiveParamGap = actualActiveTrainableParams - ...
    targetActiveTrainableParams;
stats.prunedPrunableParams = prunedPrunableParams;
stats.remainingPrunableParams = remainingPrunableParams;
stats.remainingTotalTrainableParams = remainingPrunableParams + ...
    protectedTrainableParams;
stats.actualActiveTrainableParams = actualActiveTrainableParams;
stats.actualPrunedPrunableParams = actualPrunedPrunableParams;
stats.actualPrunableSparsity = actualPrunableSparsity;
stats.activeWeightParams = activeWeightParams;
stats.activeBiasParams = activeBiasParams;
stats.prunedWeightParams = paramCounts.prunedWeightParams;
stats.prunedBiasParams = paramCounts.prunedBiasParams;
stats.numPrunedParams = actualPrunedPrunableParams;
stats.numRemainingParams = remainingPrunableParams;
stats.parameterNames = parameterNames;
stats.parameterTotal = parameterTotal;
stats.parameterPruned = parameterPruned;
stats.parameterRemaining = parameterTotal - parameterPruned;

pruningState = struct();
pruningState.masks = masks;
pruningState.parameterNames = parameterNames;
pruningState.parameterTotal = parameterTotal;
pruningState.parameterPruned = parameterPruned;
pruningState.includeBiases = pruningCfg.includeBiases;
pruningState.scope = char(scope);
pruningState.targetMode = char(pruningCfg.targetMode);
pruningState.structureMode = char(pruningCfg.structureMode);
pruningState.structuredRanking = char(pruningCfg.structuredRanking);
pruningState.structuredTargetPolicy = char(pruningCfg.structuredTargetPolicy);
pruningState.targetActiveTrainableParams = targetActiveTrainableParams;

fprintf("Pruning scope : %s\n", char(scope));
fprintf("Pruning target mode: %s\n", char(pruningCfg.targetMode));
if string(pruningCfg.targetMode) == "activeTrainableParams"
    fprintf("Target active trainable params: %d\n", ...
        targetActiveTrainableParams);
end
fprintf("Pruning target: %.2f %% effective prunable sparsity\n", ...
    100 * targetPrunableSparsity);
fprintf("Pruning actual: %.2f %% (%d/%d prunable parameters)\n", ...
    100 * actualPrunableSparsity, actualPrunedPrunableParams, ...
    totalPrunableParams);
fprintf("Active trainable params: %d/%d (weights=%d, biases=%d, protected=%d)\n", ...
    actualActiveTrainableParams, totalTrainableParams, activeWeightParams, ...
    activeBiasParams, protectedTrainableParams);
end

function pruningCfg = normalizePruningCfgForMasks(pruningCfg)
if ~isfield(pruningCfg, 'targetMode') || isempty(pruningCfg.targetMode)
    pruningCfg.targetMode = "sparsity";
end
if ~isfield(pruningCfg, 'sparsity') || isempty(pruningCfg.sparsity)
    pruningCfg.sparsity = 0;
end
if ~isfield(pruningCfg, 'targetActiveTrainableParams')
    pruningCfg.targetActiveTrainableParams = [];
end
if ~isfield(pruningCfg, 'scope') || isempty(pruningCfg.scope)
    pruningCfg.scope = "global";
end
if ~isfield(pruningCfg, 'includeBiases') || isempty(pruningCfg.includeBiases)
    pruningCfg.includeBiases = false;
end
if ~isfield(pruningCfg, 'structureMode') || isempty(pruningCfg.structureMode)
    pruningCfg.structureMode = "unstructured";
end
if ~isfield(pruningCfg, 'structuredRanking') || ...
        isempty(pruningCfg.structuredRanking)
    pruningCfg.structuredRanking = "magnitude";
end
if ~isfield(pruningCfg, 'structuredTargetPolicy') || ...
        isempty(pruningCfg.structuredTargetPolicy)
    pruningCfg.structuredTargetPolicy = "closestNotAbove";
end
if ~isfield(pruningCfg, 'hybridExactTarget') || ...
        isempty(pruningCfg.hybridExactTarget)
    pruningCfg.hybridExactTarget = false;
end

pruningCfg.targetMode = string(pruningCfg.targetMode);
pruningCfg.sparsity = double(pruningCfg.sparsity);
pruningCfg.scope = string(pruningCfg.scope);
pruningCfg.includeBiases = logical(pruningCfg.includeBiases);
pruningCfg.structureMode = string(pruningCfg.structureMode);
pruningCfg.structuredRanking = string(pruningCfg.structuredRanking);
pruningCfg.structuredTargetPolicy = string(pruningCfg.structuredTargetPolicy);
pruningCfg.hybridExactTarget = logical(pruningCfg.hybridExactTarget);
end

function [targetActiveTrainableParams, targetActivePrunableParams, ...
    prunedPrunableParams, targetPrunableSparsity] = resolveTargetCounts( ...
    pruningCfg, totalTrainableParams, totalPrunableParams, ...
    protectedTrainableParams)

if string(pruningCfg.targetMode) == "activeTrainableParams"
    targetActiveTrainableParams = double( ...
        pruningCfg.targetActiveTrainableParams);
    if isempty(targetActiveTrainableParams) || ...
            ~isscalar(targetActiveTrainableParams) || ...
            ~isfinite(targetActiveTrainableParams)
        error("cfg.pruning.targetActiveTrainableParams must be a positive integer scalar when targetMode is 'activeTrainableParams'.");
    end
    if targetActiveTrainableParams < protectedTrainableParams
        error("Requested active trainable parameter target is smaller than the protected parameter count. Either increase the target or enable bias pruning / include more parameters in the prunable set.");
    end
    if targetActiveTrainableParams > totalTrainableParams
        error("Requested active trainable parameter target exceeds total trainable parameter count.");
    end

    targetActivePrunableParams = targetActiveTrainableParams - ...
        protectedTrainableParams;
    prunedPrunableParams = totalPrunableParams - targetActivePrunableParams;
else
    prunedPrunableParams = floor(double(pruningCfg.sparsity) * ...
        totalPrunableParams);
    targetActivePrunableParams = totalPrunableParams - ...
        prunedPrunableParams;
    targetActiveTrainableParams = targetActivePrunableParams + ...
        protectedTrainableParams;
end

prunedPrunableParams = max(0, min(totalPrunableParams, ...
    round(prunedPrunableParams)));
targetActivePrunableParams = totalPrunableParams - prunedPrunableParams;
targetPrunableSparsity = prunedPrunableParams / max(totalPrunableParams, 1);
end

function [masks, parameterNames, parameterTotal, parameterPruned] = ...
    createGlobalMasks(masks, candidates, numToPrune)
allMagnitudes = vertcat(candidates.magnitudes);
pruneFlags = false(numel(allMagnitudes), 1);
if numToPrune > 0
    [~, order] = sort(allMagnitudes, "ascend");
    pruneFlags(order(1:numToPrune)) = true;
end

offset = 0;
parameterNames = strings(numel(candidates), 1);
parameterTotal = zeros(numel(candidates), 1);
parameterPruned = zeros(numel(candidates), 1);

for i = 1:numel(candidates)
    idx = offset + (1:candidates(i).numel);
    keepMask = reshape(~pruneFlags(idx), candidates(i).size);
    masks{candidates(i).row} = keepMask;

    parameterNames(i) = candidates(i).name;
    parameterTotal(i) = candidates(i).numel;
    parameterPruned(i) = nnz(~keepMask);
    offset = offset + candidates(i).numel;
end
end

function [masks, parameterNames, parameterTotal, parameterPruned] = ...
    createLayerwiseMasks(masks, candidates, totalNumToPrune)
parameterNames = strings(numel(candidates), 1);
parameterTotal = zeros(numel(candidates), 1);
parameterPruned = zeros(numel(candidates), 1);
if isempty(candidates)
    return;
end

counts = [candidates.numel];
prunePerTensor = allocateLayerwisePrunedCounts(totalNumToPrune, counts, ...
    string({candidates.name}));

for i = 1:numel(candidates)
    magnitudes = candidates(i).magnitudes;
    numToPrune = prunePerTensor(i);
    pruneFlags = false(numel(magnitudes), 1);
    if numToPrune > 0
        [~, order] = sort(magnitudes, "ascend");
        pruneFlags(order(1:numToPrune)) = true;
    end

    keepMask = reshape(~pruneFlags, candidates(i).size);
    masks{candidates(i).row} = keepMask;

    parameterNames(i) = candidates(i).name;
    parameterTotal(i) = candidates(i).numel;
    parameterPruned(i) = nnz(~keepMask);
end
end

function prunePerTensor = allocateLayerwisePrunedCounts(totalNumToPrune, ...
    counts, names)
counts = double(counts(:).');
totalNumToPrune = round(double(totalNumToPrune));
if totalNumToPrune <= 0
    prunePerTensor = zeros(size(counts));
    return;
end

raw = totalNumToPrune * counts / sum(counts);
prunePerTensor = floor(raw);
fractional = raw - prunePerTensor;
remaining = totalNumToPrune - sum(prunePerTensor);

[~, order] = sortrows(table(-fractional(:), -counts(:), names(:)), ...
    [1 2 3]);
while remaining > 0
    madeProgress = false;
    for k = 1:numel(order)
        idx = order(k);
        if prunePerTensor(idx) < counts(idx)
            prunePerTensor(idx) = prunePerTensor(idx) + 1;
            remaining = remaining - 1;
            madeProgress = true;
            if remaining == 0
                break;
            end
        end
    end
    if ~madeProgress
        break;
    end
end
end

function tf = isPrunableParameter(parameterName, includeBiases)
name = lower(char(string(parameterName)));
tf = strcmp(name, "weights") || (includeBiases && strcmp(name, "bias"));
end

function name = learnableName(learnables, row)
name = string(learnables.Layer(row)) + "/" + string(learnables.Parameter(row));
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    data = extractdata(value);
else
    data = value;
end
data = gather(data);
end
