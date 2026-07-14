function counts = summarizeTrainableParameters(net, includeBiases, masks)
% summarizeTrainableParameters - Count active and pruned trainable parameters.
%
% Counts weights, biases, protected parameters, and prunable parameters from a
% dlnetwork. Optional masks let dense and pruned networks share the same
% accounting path.

if nargin < 2 || isempty(includeBiases)
    includeBiases = false;
end
includeBiases = logical(includeBiases);

learnables = net.Learnables;
if nargin < 3 || isempty(masks)
    masks = cell(height(learnables), 1);
    for row = 1:height(learnables)
        masks{row} = true(size(learnableToNumeric(learnables.Value{row})));
    end
end

counts = struct();
counts.totalTrainableParams = 0;
counts.totalPrunableParams = 0;
counts.protectedTrainableParams = 0;
counts.totalWeightParams = 0;
counts.totalBiasParams = 0;
counts.activeWeightParams = 0;
counts.activeBiasParams = 0;
counts.activeOtherParams = 0;
counts.prunedWeightParams = 0;
counts.prunedBiasParams = 0;
counts.prunedOtherParams = 0;
counts.actualActiveTrainableParams = 0;
counts.actualPrunedPrunableParams = 0;
counts.remainingPrunableParams = 0;
counts.remainingTotalTrainableParams = 0;
counts.actualPrunableSparsity = 0;

for row = 1:height(learnables)
    data = learnableToNumeric(learnables.Value{row});
    totalCount = numel(data);
    if row <= numel(masks) && ~isempty(masks{row})
        keepMask = logical(masks{row});
    else
        keepMask = true(size(data));
    end
    activeCount = nnz(keepMask);
    prunedCount = totalCount - activeCount;

    parameterName = lower(string(learnables.Parameter(row)));
    isWeight = parameterName == "weights";
    isBias = parameterName == "bias";
    isPrunable = isWeight || (includeBiases && isBias);

    counts.totalTrainableParams = counts.totalTrainableParams + totalCount;
    counts.actualActiveTrainableParams = ...
        counts.actualActiveTrainableParams + activeCount;

    if isPrunable
        counts.totalPrunableParams = counts.totalPrunableParams + totalCount;
        counts.actualPrunedPrunableParams = ...
            counts.actualPrunedPrunableParams + prunedCount;
    else
        counts.protectedTrainableParams = ...
            counts.protectedTrainableParams + totalCount;
    end

    if isWeight
        counts.totalWeightParams = counts.totalWeightParams + totalCount;
        counts.activeWeightParams = counts.activeWeightParams + activeCount;
        counts.prunedWeightParams = counts.prunedWeightParams + prunedCount;
    elseif isBias
        counts.totalBiasParams = counts.totalBiasParams + totalCount;
        counts.activeBiasParams = counts.activeBiasParams + activeCount;
        counts.prunedBiasParams = counts.prunedBiasParams + prunedCount;
    else
        counts.activeOtherParams = counts.activeOtherParams + activeCount;
        counts.prunedOtherParams = counts.prunedOtherParams + prunedCount;
    end
end

counts.remainingPrunableParams = counts.totalPrunableParams - ...
    counts.actualPrunedPrunableParams;
counts.remainingTotalTrainableParams = counts.actualActiveTrainableParams;
counts.actualPrunableSparsity = counts.actualPrunedPrunableParams / ...
    max(counts.totalPrunableParams, 1);
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    data = extractdata(value);
else
    data = value;
end
data = gather(data);
end
