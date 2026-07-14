function yhat = predictPhaseNorm(netDPD, inputMtxN, normStats, r_sel)
% predictPhaseNorm - Reconstruct complex PNNN output from normalized NN predictions.
%
% This helper runs the trained phase-normalized network, denormalizes the
% two-channel prediction, and rotates it back to the complex signal domain
% used by the offline evaluation flow.
%
% Inputs:
%   netDPD - Trained PNNN network.
%   inputMtxN - Normalized input feature matrix.
%   normStats - Struct with output normalization fields muY and sigmaY.
%   r_sel - Phase-rotation vector aligned with the selected samples.
%
% Outputs:
%   yhat - Reconstructed complex model output.
%
% Notes:
%   This function does not decide the physical X/Y semantics; it only
%   reconstructs the complex output from the normalized phase-domain prediction.

predN = predict(netDPD, inputMtxN);
pred = predN .* normStats.sigmaY + normStats.muY;
y_rot = pred(:,1) + 1j*pred(:,2);
yhat = conj(r_sel(:)) .* y_rot(:);
end
