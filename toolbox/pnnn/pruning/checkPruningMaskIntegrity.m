function checkPruningMaskIntegrity(net, masks)
% checkPruningMaskIntegrity - Require every masked PNNN value to remain zero.

learnables = net.Learnables;
for row = 1:height(learnables)
    value = learnables.Value{row};
    if isa(value, 'dlarray')
        value = extractdata(value);
    end
    value = gather(value);
    if any(value(~logical(masks{row})) ~= 0, 'all')
        error('checkPruningMaskIntegrity:MaskViolation', ...
            'A pruned PNNN parameter became nonzero.');
    end
end
end
