function counts = summarizeTrainableParameters(net, masks)
% summarizeTrainableParameters - Count total and active PNNN weights/biases.

learnables = net.Learnables;
if nargin < 2
    masks = cell(height(learnables), 1);
    for row = 1:height(learnables)
        masks{row} = true(size(learnableToNumeric(learnables.Value{row})));
    end
end

counts = struct('totalWeightParams', 0, 'totalBiasParams', 0, ...
    'activeWeightParams', 0, 'activeBiasParams', 0);

for row = 1:height(learnables)
    total = numel(learnableToNumeric(learnables.Value{row}));
    active = nnz(masks{row});
    name = lower(string(learnables.Parameter(row)));
    if name == "weights"
        counts.totalWeightParams = counts.totalWeightParams + total;
        counts.activeWeightParams = counts.activeWeightParams + active;
    elseif name == "bias"
        counts.totalBiasParams = counts.totalBiasParams + total;
        counts.activeBiasParams = counts.activeBiasParams + active;
    end
end
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    value = extractdata(value);
end
data = gather(value);
end
