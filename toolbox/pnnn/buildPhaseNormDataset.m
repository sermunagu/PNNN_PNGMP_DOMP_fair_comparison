function [features, neuralTargets, phaseRotation] = ...
    buildPhaseNormDataset(x, y, M, orders, featMode, options)
% buildPhaseNormDataset - Build the complete phase-normalized PNNN dataset.
% Memory taps use periodic extension: x(n), x(n-1), ..., x(n-M).

if nargin < 5 || isempty(featMode)
    featMode = "full";
end

x = x(:).';
y = y(:).';
orders = orders(:).';
featMode = string(featMode);
sampleCount = numel(x);

switch featMode
    case "full"
        featureCount = 2*(M+1) + (M+1)*numel(orders);
    case "pruned"
        featureCount = 1 + 2*M + M*numel(orders);
    otherwise
        error("featMode debe ser 'full' o 'pruned'.");
end

features = zeros(featureCount, sampleCount);
phaseRotation = complex(zeros(1, sampleCount));

for sample = 1:sampleCount
    if abs(x(sample)) == 0
        rotation = 1 + 0j;
    else
        rotation = conj(x(sample)) / abs(x(sample));
    end

    taps = mod((sample:-1:sample-M) - 1, sampleCount) + 1;
    tappedInput = x(taps);
    rotatedTaps = rotation * tappedInput;

    switch featMode
        case "full"
            realTaps = real(rotatedTaps);
            imaginaryTaps = imag(rotatedTaps);
            amplitudes = abs(tappedInput);

            envelopes = zeros(1, (M+1)*numel(orders));
            for orderIndex = 1:numel(orders)
                columns = (orderIndex-1)*(M+1) + (1:M+1);
                envelopes(columns) = amplitudes.^orders(orderIndex);
            end
            featureVector = [realTaps, imaginaryTaps, envelopes];

        case "pruned"
            realMainTap = real(rotatedTaps(1));
            realMemoryTaps = real(rotatedTaps(2:end));
            imaginaryMemoryTaps = imag(rotatedTaps(2:end));
            memoryAmplitudes = abs(tappedInput(2:end));

            envelopes = zeros(1, M*numel(orders));
            for orderIndex = 1:numel(orders)
                columns = (orderIndex-1)*M + (1:M);
                envelopes(columns) = ...
                    memoryAmplitudes.^orders(orderIndex);
            end
            featureVector = [realMainTap, realMemoryTaps, ...
                imaginaryMemoryTaps, envelopes];
    end

    features(:, sample) = featureVector(:);
    phaseRotation(sample) = rotation;
end

% Explicit experiment only: preserve all 84 positions and the historical default.
if nargin >= 6 && isfield(options, 'zeroCurrentImaginaryTap') && ...
        options.zeroCurrentImaginaryTap && featMode == "full"
    features(M+2,:) = 0;
end
neuralTargets = zeros(2, sampleCount);
for sample = 1:sampleCount
    rotatedTarget = phaseRotation(sample) * y(sample);
    neuralTargets(:, sample) = [real(rotatedTarget); imag(rotatedTarget)];
end
end
