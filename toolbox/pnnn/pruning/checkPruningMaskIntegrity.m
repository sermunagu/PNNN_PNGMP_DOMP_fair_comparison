function [integrity, stats] = checkPruningMaskIntegrity(net, pruningState, stats, stage)
% checkPruningMaskIntegrity - Verify that pruned PNNN weights stay zero.
%
% This function inspects masked learnables, counts any nonzero values at
% pruned positions, and updates pruning metadata after pruning/fine-tuning.
%
% Inputs:
%   net - dlnetwork to inspect.
%   pruningState - Struct containing the masks aligned with net.Learnables.
%   stats - Existing pruning statistics struct.
%   stage - Text label for the integrity check stage.
%
% Outputs:
%   integrity - Struct with violation count, maximum magnitude, and status.
%   stats - Updated pruning statistics struct.

if nargin < 4
    stage = "";
end

if ~isa(net, 'dlnetwork')
    error("checkPruningMaskIntegrity requiere un objeto dlnetwork.");
end
if ~isfield(pruningState, 'masks')
    error("pruningState debe contener el campo masks.");
end

tolerance = 0;
learnables = net.Learnables;
masks = pruningState.masks;
if numel(masks) ~= height(learnables)
    error("El numero de mascaras no coincide con net.Learnables.");
end

violationCount = 0;
violationMaxAbs = 0;

for i = 1:height(learnables)
    mask = masks{i};
    if isempty(mask)
        continue;
    end

    data = learnableToNumeric(learnables.Value{i});
    if ~isequal(size(mask), size(data))
        error("La mascara %d no coincide con el tamano del parametro learnable.", i);
    end

    prunedValues = abs(data(~logical(mask)));
    if isempty(prunedValues)
        continue;
    end

    violationCount = violationCount + nnz(prunedValues > tolerance);
    violationMaxAbs = max(violationMaxAbs, max(prunedValues(:)));
end

integrity = struct();
integrity.stage = char(string(stage));
integrity.violationCount = violationCount;
integrity.violationMaxAbs = violationMaxAbs;
integrity.tolerance = tolerance;
integrity.ok = (violationCount == 0);

if nargin >= 3 && ~isempty(stats)
    stats.maskViolationCount = integrity.violationCount;
    stats.maskViolationMaxAbs = integrity.violationMaxAbs;
    stats.maskIntegrityOk = integrity.ok;
    stats.maskIntegrityStage = string(stage);
    stats.maskIntegrityChecks(end+1) = integrity;
end
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    data = extractdata(value);
else
    data = value;
end
data = gather(data);
end
