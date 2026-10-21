function counts = summarizeTrainableParameters(network, masks)
% summarizeTrainableParameters - Count total and active PNNN parameters.

learnables = network.Learnables;
counts = struct('totalWeightParams', 0, 'totalBiasParams', 0, ...
    'activeWeightParams', 0, 'activeBiasParams', 0, ...
    'maxAbsRealParameter', 0);

for row = 1:height(learnables)
    data = gather(extractdata(learnables.Value{row}));
    if nargin < 2
        mask = true(size(data));
    else
        mask = logical(masks{row});
    end
    if ~isequal(size(mask), size(data))
        error('summarizeTrainableParameters:MaskSizeMismatch', ...
            'Each learnable mask must match its parameter array.');
    end
    total = numel(data);
    active = nnz(mask);
    
    if active > 0
        counts.maxAbsRealParameter = max(counts.maxAbsRealParameter, max(abs(data(mask))));
    end
    name = lower(string(learnables.Parameter(row)));
    if name == "weights"
        counts.totalWeightParams = counts.totalWeightParams + total;
        counts.activeWeightParams = counts.activeWeightParams + active;
    elseif name == "bias"
        counts.totalBiasParams = counts.totalBiasParams + total;
        counts.activeBiasParams = counts.activeBiasParams + active;
    end
end
counts.maxAbsRealParameter = double(counts.maxAbsRealParameter);
end
