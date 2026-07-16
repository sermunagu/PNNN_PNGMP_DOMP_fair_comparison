function description = describeAdditionalOperations(row, activationName)
% describeAdditionalOperations - Human-readable special operations per sample.

if nargin < 2
    activationName = "activation";
end
parts = strings(0, 1);

if row.NumAbsPerSample > 0
    if row.NumSqrtPerSample == row.NumAbsPerSample
        parts(end+1) = sprintf('%d magnitude (including %d sqrt)', ...
            row.NumAbsPerSample, row.NumSqrtPerSample);
    else
        parts(end+1) = sprintf('%d magnitude', row.NumAbsPerSample);
    end
end
if row.NumSqrtPerSample > 0 && row.NumSqrtPerSample ~= row.NumAbsPerSample
    parts(end+1) = sprintf('%d sqrt', row.NumSqrtPerSample);
end
if row.NumRealDivisionsPerSample > 0
    parts(end+1) = sprintf('%d division', row.NumRealDivisionsPerSample);
end
if row.NumActivationEvaluationsPerSample > 0
    parts(end+1) = sprintf('%d %s', ...
        row.NumActivationEvaluationsPerSample, activationName);
end
if row.NumExpWorstCasePerSample > 0
    parts(end+1) = sprintf('up to %d exp', row.NumExpWorstCasePerSample);
end
if row.PhaseNormalizationIncluded
    parts(end+1) = ...
        "phase normalization/restoration arithmetic in FLOPs";
end

if isempty(parts)
    description = "none";
else
    description = strjoin(parts, ', ');
end
end
