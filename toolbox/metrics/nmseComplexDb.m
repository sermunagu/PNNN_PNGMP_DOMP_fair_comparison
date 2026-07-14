function value = nmseComplexDb(target, prediction)
% nmseComplexDb - Compute the shared complex time-domain NMSE in dB.
% Every polynomial and neural model in the fair comparison uses this exact
% energy-ratio definition on the same explicitly supplied sample domain.

target = double(target(:));
prediction = double(prediction(:));
if isempty(target) || numel(prediction) ~= numel(target) || ...
        any(~isfinite(target)) || any(~isfinite(prediction))
    error('nmseComplexDb:InvalidSignals', ...
        'target and prediction must be finite vectors of equal length.');
end
target_energy = sum(abs(target).^2);
if target_energy <= 0 || ~isfinite(target_energy)
    error('nmseComplexDb:ZeroTargetEnergy', ...
        'The target must have positive finite energy.');
end
error_energy = sum(abs(target - prediction).^2);
value = 10*log10(max(error_energy, realmin)/target_energy);
end
