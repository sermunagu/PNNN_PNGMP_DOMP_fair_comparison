function features = buildPNIQFeatures( ...
    x, rows, rotation, manager, support, descriptors)
% Build the real I/Q representation of phase-normalized GMP regressors.

regressorsI = zeros(numel(rows), numel(support));
regressorsQ = zeros(numel(rows), numel(support));
signalLength = numel(x);

for localIndex = 1:numel(support)
    descriptor = descriptors(support(localIndex));
    if descriptor.canonicalGMP
        carrierRows = mod( ...
            rows - descriptor.carrierLag - 1, signalLength) + 1;
        normalizedCarrier = rotation .* x(carrierRows);
        envelope = ones(numel(rows), 1);
        for termIndex = 1:numel(descriptor.envelopeLags)
            envelopeRows = mod(rows - ...
                descriptor.envelopeLags(termIndex) - 1, signalLength) + 1;
            envelope = envelope .* abs(x(envelopeRows)).^ ...
                descriptor.envelopePowers(termIndex);
        end
        regressorsI(:, localIndex) = real(normalizedCarrier) .* envelope;
        regressorsQ(:, localIndex) = imag(normalizedCarrier) .* envelope;
        if descriptor.QColumnStructurallyZero
            regressorsQ(:, localIndex) = 0;
        end
    else
        % Only auxiliaries need the unfactorized complex regressor.
        regressor = buildGMPRegressorRows(x, rows, manager, support(localIndex));
        phaseNormalized = rotation .* regressor;
        regressorsI(:, localIndex) = real(phaseNormalized);
        regressorsQ(:, localIndex) = imag(phaseNormalized);
    end
end
features = [regressorsI, regressorsQ];
end
