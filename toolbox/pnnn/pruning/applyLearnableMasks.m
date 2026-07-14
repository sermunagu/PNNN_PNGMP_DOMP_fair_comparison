function net = applyLearnableMasks(net, masks)
% applyLearnableMasks - Apply binary masks to dlnetwork learnables.
%
% This function multiplies each learnable value by its mask and writes the
% masked learnables back to the network during pruning and fine-tuning.
%
% Inputs:
%   net - dlnetwork whose learnables must be masked.
%   masks - Cell array aligned with net.Learnables rows.
%
% Outputs:
%   net - Network with masked learnable values.

if ~isa(net, 'dlnetwork')
    error("applyLearnableMasks requiere un objeto dlnetwork.");
end

learnables = net.Learnables;
if numel(masks) ~= height(learnables)
    error("El numero de mascaras no coincide con net.Learnables.");
end

for i = 1:height(learnables)
    if ~isempty(masks{i})
        learnables.Value{i} = applyMaskToValue(learnables.Value{i}, masks{i});
    end
end

net = assignLearnables(net, learnables);
end

function value = applyMaskToValue(value, mask)
mask = logical(mask);

if isa(value, 'dlarray')
    data = extractdata(value);
    maskValue = double(mask);
    if isa(data, 'gpuArray')
        maskValue = gpuArray(maskValue);
    end
    maskValue = maskValue .* ones(size(data), "like", data);
    value = value .* maskValue;
else
    maskValue = double(mask);
    if isa(value, 'gpuArray')
        maskValue = gpuArray(maskValue);
    end
    maskValue = maskValue .* ones(size(value), "like", value);
    value = value .* maskValue;
end
end

function net = assignLearnables(net, learnables)
try
    net.Learnables = learnables;
catch ME
    if exist("setLearnableParameterValue", "file") ~= 2
        rethrow(ME);
    end

    for i = 1:height(learnables)
        net = setLearnableParameterValue(net, ...
            char(string(learnables.Layer(i))), ...
            char(string(learnables.Parameter(i))), ...
            learnables.Value{i});
    end
end
end
