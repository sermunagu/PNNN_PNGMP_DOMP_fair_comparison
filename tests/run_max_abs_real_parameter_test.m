% Verify peak-normalized maxima for linear families and native PNNN.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));

%% Complex GMP definition and prediction-coefficient transformation
selectedRegressors = [ ...
    1+2i, 2-1i; 2-1i, -1+3i; -1+1i, 3+2i; ...
    3-2i, 1+1i; 2+3i, -2+1i; -2-1i, 1-3i];
identificationTarget = [2+1i; -1+4i; 3-2i; -2-3i; 1+5i; 4-1i];
regressorPeaks = max(abs(selectedRegressors), [], 1).';
unitPeakRegressors = selectedRegressors ./ regressorPeaks.';
assert(max(abs(max(abs(unitPeakRegressors), [], 1) - 1)) < 1e-14);
outputPeak = max(abs(identificationTarget));
unitPeakIdentificationTarget = identificationTarget / outputPeak;
tolerance = max(size(unitPeakRegressors)) * ...
    eps(norm(unitPeakRegressors, 2));
expectedCoefficients = lsqminnorm( ...
    unitPeakRegressors, unitPeakIdentificationTarget, tolerance);
expectedMaximum = max(abs(expectedCoefficients(:)));
implementedMaximum = max(abs(expectedCoefficients(:)));
assert(abs(implementedMaximum - expectedMaximum) < 1e-14);
assert(abs(abs(3 + 4i) - 5) < 1e-14);
assert(max([abs(real(3 + 4i)); abs(imag(3 + 4i))]) ~= abs(3 + 4i));

predictionCoefficients = ...
    (expectedCoefficients * outputPeak) ./ regressorPeaks;
predictionFromOriginal = selectedRegressors * predictionCoefficients;
predictionFromNormalized = ...
    outputPeak * unitPeakRegressors * expectedCoefficients;
assert(max(abs(predictionFromOriginal - predictionFromNormalized)) < 1e-12);

scaledTarget = 2.7 * identificationTarget;
scaledOutputPeak = max(abs(scaledTarget));
scaledCoefficients = lsqminnorm(unitPeakRegressors, ...
    scaledTarget/scaledOutputPeak, tolerance);
assert(max(abs(scaledCoefficients - expectedCoefficients)) < 1e-12);
columnScales = [2.5 0.4];
scaledRegressors = selectedRegressors .* columnScales;
scaledRegressorPeaks = max(abs(scaledRegressors), [], 1).';
assert(max(abs(scaledRegressors ./ scaledRegressorPeaks.' - ...
    unitPeakRegressors), [], 'all') < 1e-14);

%% PN-IQ definition and prediction-coefficient transformation
selectedFeatures = [1 2; -2 1; 3 -1; 1 -3; -1 2; 2 1];
angles = (0:size(selectedFeatures, 1)-1).' * pi/11;
identificationRotation = exp(1i*angles);
phaseNormalizedIdentificationTarget = ...
    identificationRotation .* identificationTarget;
outputPeakPNIQ = max(abs(phaseNormalizedIdentificationTarget));
unitPeakPhaseNormalizedIdentificationTarget = ...
    phaseNormalizedIdentificationTarget / outputPeakPNIQ;
featurePeaks = max(abs(selectedFeatures), [], 1).';
unitPeakFeatures = selectedFeatures ./ featurePeaks.';
assert(max(abs(max(abs(unitPeakFeatures), [], 1) - 1)) < 1e-14);
tolerancePNIQ = max(size(unitPeakFeatures)) * ...
    eps(norm(unitPeakFeatures, 2));
expectedI = lsqminnorm(unitPeakFeatures, ...
    real(unitPeakPhaseNormalizedIdentificationTarget), tolerancePNIQ);
expectedQ = lsqminnorm(unitPeakFeatures, ...
    imag(unitPeakPhaseNormalizedIdentificationTarget), tolerancePNIQ);
expectedPNIQMaximum = max(abs([expectedI(:); expectedQ(:)]));
implementedPNIQMaximum = max(abs([expectedI(:); expectedQ(:)]));
assert(abs(implementedPNIQMaximum - expectedPNIQMaximum) < 1e-14);

predictionCoefficientsI = ...
    (expectedI * outputPeakPNIQ) ./ featurePeaks;
predictionCoefficientsQ = ...
    (expectedQ * outputPeakPNIQ) ./ featurePeaks;
assert(max(abs(selectedFeatures*predictionCoefficientsI - ...
    outputPeakPNIQ*unitPeakFeatures*expectedI)) < 1e-12);
assert(max(abs(selectedFeatures*predictionCoefficientsQ - ...
    outputPeakPNIQ*unitPeakFeatures*expectedQ)) < 1e-12);

scaledPhaseTarget = 2.7 * phaseNormalizedIdentificationTarget;
scaledPeakPNIQ = max(abs(scaledPhaseTarget));
scaledI = lsqminnorm(unitPeakFeatures, ...
    real(scaledPhaseTarget/scaledPeakPNIQ), tolerancePNIQ);
scaledQ = lsqminnorm(unitPeakFeatures, ...
    imag(scaledPhaseTarget/scaledPeakPNIQ), tolerancePNIQ);
assert(max(abs([scaledI-expectedI; scaledQ-expectedQ])) < 1e-12);

%% Fixed-Ridge definitions
lambda = 1e-3;
gram = unitPeakRegressors' * unitPeakRegressors;
rhs = unitPeakRegressors' * unitPeakIdentificationTarget;
expectedRidgeCoefficients = ...
    (gram + lambda*eye(size(gram))) \ rhs;
assert(max(abs(expectedRidgeCoefficients(:))) > 0);

gramPNIQ = unitPeakFeatures' * unitPeakFeatures;
rhsI = unitPeakFeatures' * real(unitPeakPhaseNormalizedIdentificationTarget);
rhsQ = unitPeakFeatures' * imag(unitPeakPhaseNormalizedIdentificationTarget);
expectedRidgeI = (gramPNIQ + lambda*eye(size(gramPNIQ))) \ rhsI;
expectedRidgeQ = (gramPNIQ + lambda*eye(size(gramPNIQ))) \ rhsQ;
assert(max(abs([expectedRidgeI(:); expectedRidgeQ(:)])) > 0);

%% PNNN active-mask definition remains unchanged
layers = [ ...
    featureInputLayer(2, 'Normalization', 'none', 'Name', 'input')
    fullyConnectedLayer(1, 'Name', 'output')];
network = dlnetwork(layers);
learnables = network.Learnables;
weightRow = lower(string(learnables.Parameter)) == "weights";
biasRow = lower(string(learnables.Parameter)) == "bias";
assert(nnz(weightRow) == 1 && nnz(biasRow) == 1);
learnables.Value{weightRow} = dlarray([100 3]);
learnables.Value{biasRow} = dlarray(4);
network.Learnables = learnables;
masks = cell(height(learnables), 1);
for row = 1:height(learnables)
    masks{row} = true(size(extractdata(learnables.Value{row})));
end
weightMask = masks{weightRow};
weightMask(1) = false;
masks{weightRow} = weightMask;
counts = summarizeTrainableParameters(network, masks);
assert(counts.activeWeightParams == 1);
assert(counts.activeBiasParams == 1);
assert(counts.maxAbsRealParameter == 4);

%% Source and presentation contracts
complexSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_complex_gmp_domp.m'));
pnSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_pniq_gmp.m'));
fixedSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'run_fixed_ridge_sweep.m'));
sweepSource = fileread(fullfile(projectRoot, 'run_parameter_sweep.m'));
assert(contains(complexSource, 'unitPeakRegressionCoefficientPaths'));
assert(contains(pnSource, 'unitPeakRegressionCoefficientPathsI'));
assert(contains(pnSource, 'unitPeakRegressionCoefficientPathsQ'));
assert(contains(fixedSource, 'complexUnitPeakRegressionCoefficientPaths'));
assert(contains(fixedSource, 'pniqUnitPeakRegressionCoefficientPathsI'));
assert(contains(fixedSource, 'pniqUnitPeakRegressionCoefficientPathsQ'));
assert(isscalar(strfind(complexSource, 'lsqminnorm')));
assert(numel(strfind(pnSource, 'lsqminnorm')) == 2);
assert(~contains(complexSource, 'diagLoad'));
assert(~contains(pnSource, 'diagLoad'));
assert(~contains(fixedSource, 'diagLoad'));

cfg = getFairDOMPComparisonConfig(projectRoot);
assert(cfg.sweep.coefficientRangeDefinition == ...
    "unit_peak_output_per_column_peak_regressors_v4");
assert(all(strlength([cfg.names.complexGMPDOMP; ...
    cfg.names.pniqGMP; cfg.names.pnnn]) > 0));
assert(contains(sweepSource, ...
    "'metricLabel', 'Maximum absolute coefficient/parameter'"));
assert(contains(sweepSource, "'names', cfg.names"));

fprintf('MAX ABS REAL PARAMETER TEST: PASS\n');
