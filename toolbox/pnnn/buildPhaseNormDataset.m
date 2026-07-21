function [features, neuralTargets, phaseRotation] = ...
    buildPhaseNormDataset(x, y, M, orders, featMode)
% buildPhaseNormDataset - Build the complete phase-normalized PNNN dataset.
% Memory taps use periodic extension: x(n), x(n-1), ..., x(n-M).

%% 1. Validate and align the modeled input and output signals
if nargin < 5 || isempty(featMode)
    featMode = "full";
end

x = x(:).';
y = y(:).';
if numel(x) ~= numel(y)
    error('x e y deben tener la misma longitud.');
end
if any(~isfinite(x)) || any(~isfinite(y))
    error('x o y contienen NaN/Inf.');
end

%% 2. Rotate the input and build the periodic nonlinear features
[features, phaseRotation] = buildPhaseNormInput(x, M, orders, featMode);

%% 3. Rotate the complex target into two real output channels
sampleCount = numel(x);
neuralTargets = zeros(2, sampleCount);
for sample = 1:sampleCount
    rotatedTarget = phaseRotation(sample) * y(sample);
    neuralTargets(:, sample) = [real(rotatedTarget); imag(rotatedTarget)];
end
end

function [features, phaseRotation] = ...
    buildPhaseNormInput(x, M, orders, featMode)
% Build the feature matrix without changing the established feature order.

x = x(:).';
validateInputSignal(x, M);
orders = validateOrders(orders);
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

    taps = periodicTapIndices(sample, M, sampleCount);
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
end

function taps = periodicTapIndices(sample, M, sampleCount)
taps = mod((sample:-1:sample-M) - 1, sampleCount) + 1;
end

function validateInputSignal(x, M)
if isempty(x)
    error('x no puede estar vacío.');
end
if any(~isfinite(x))
    error('x contiene NaN/Inf.');
end
if ~isscalar(M) || M < 0 || M ~= floor(M)
    error('M debe ser un entero no negativo.');
end
if numel(x) <= M
    error('La señal debe tener más muestras que M. N=%d, M=%d.', ...
        numel(x), M);
end
end

function orders = validateOrders(orders)
orders = orders(:).';
if isempty(orders)
    error('orders no puede estar vacío.');
end
if any(~isfinite(orders)) || any(orders < 1) || ...
        any(orders ~= floor(orders))
    error('orders debe contener enteros positivos finitos.');
end
end
