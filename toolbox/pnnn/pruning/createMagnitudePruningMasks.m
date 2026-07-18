function masks = createMagnitudePruningMasks(net, targetActiveParameters)
% createMagnitudePruningMasks - Prune the globally smallest PNNN weights.
% Biases are protected; every target is an exact active-parameter count.

learnables = net.Learnables;
masks = cell(height(learnables), 1);
weightRows = zeros(0, 1);
weightSizes = cell(0, 1);
magnitudes = cell(0, 1);
totalParameters = 0;
totalWeights = 0;

for row = 1:height(learnables)
    value = learnableToNumeric(learnables.Value{row});
    masks{row} = true(size(value));
    totalParameters = totalParameters + numel(value);
    if lower(string(learnables.Parameter(row))) == "weights"
        weightRows(end+1, 1) = row; %#ok<AGROW>
        weightSizes{end+1, 1} = size(value); %#ok<AGROW>
        magnitudes{end+1, 1} = abs(value(:)); %#ok<AGROW>
        totalWeights = totalWeights + numel(value);
    end
end

protectedParameters = totalParameters - totalWeights;
if targetActiveParameters < protectedParameters || ...
        targetActiveParameters > totalParameters
    error('createMagnitudePruningMasks:InvalidTarget', ...
        'The active-parameter target is incompatible with protected biases.');
end

pruneCount = totalParameters - targetActiveParameters;
prune = false(totalWeights, 1);
if pruneCount > 0
    allMagnitudes = vertcat(magnitudes{:});
    [~, order] = sort(allMagnitudes, 'ascend');
    prune(order(1:pruneCount)) = true;
end

offset = 0;
for index = 1:numel(weightRows)
    count = numel(magnitudes{index});
    local = offset + (1:count);
    masks{weightRows(index)} = reshape(~prune(local), weightSizes{index});
    offset = offset + count;
end
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    value = extractdata(value);
end
data = gather(value);
end
