function yhat = predictPNNN(network, normalizedFeatures, normalization, phaseRotation)
% Predict, undo the stored output z-score, and restore the input phase.

predN = predict(network, normalizedFeatures);
pred = predN .* normalization.sigmaY + normalization.muY;
y_rot = pred(:,1) + 1j*pred(:,2);
yhat = conj(phaseRotation(:)) .* y_rot(:);
end
