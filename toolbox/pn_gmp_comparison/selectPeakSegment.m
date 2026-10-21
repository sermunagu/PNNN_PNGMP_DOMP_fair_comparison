function indices = selectPeakSegment(~, target, fraction)
% Historical target-peak segment selector; the unused input keeps the API.

[~, peakIndex] = max(abs(target));
firstIndex = floor(length(target)*fraction*floor((peakIndex/length(target))/fraction))+1;
lastIndex = ceil(length(target)*fraction*ceil((peakIndex/length(target))/fraction))-1;

if(lastIndex>length(target))
    lastIndex = length(target);
    firstIndex = length(target)-floor(fraction*length(target));
end

indices = firstIndex:lastIndex;
