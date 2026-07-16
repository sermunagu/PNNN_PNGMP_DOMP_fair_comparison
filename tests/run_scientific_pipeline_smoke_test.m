% Exercise the shared split, GMP regressors, PN-IQ rotation, and PNNN features.
% The deterministic fixture states the core scientific contracts directly.
% It performs no neural training and uses only small linear fits.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
for folder = ["domp","metrics","pn_gmp_comparison","pnnn","splits"]
    addpath(fullfile(projectRoot, 'toolbox', folder));
end

rng(712, 'twister');
n = 640;
x = (randn(n, 1) + 1j*randn(n, 1))/sqrt(2);
x(1:97:end) = 0;
y = 0.82*x + 0.13*circshift(x, 1).*abs(circshift(x, 2)).^2 + ...
    0.04*conj(circshift(x, 3));
cfg = getFairDOMPComparisonConfig(projectRoot);

%% The current split remains deterministic and unchanged
split = buildCommonComparisonSplit(x, y, cfg);
repeat = buildCommonComparisonSplit(x, y, cfg);
assert(isequaln(split, repeat));
assert(abs(numel(split.identificationIndices) - ...
    floor(cfg.identificationFraction*n)) <= 2);
assert(isempty(intersect(split.internalTrainIndices, ...
    split.internalValidationIndices)));
assert(isequal(sort([split.internalTrainIndices; ...
    split.internalValidationIndices]), sort(split.identificationIndices)));
assert(isequal(split.fullSignalIndices, (1:n).'));

%% Phase normalization uses the current input sample and preserves zeros
rows = split.identificationIndices;
rotation = computePhaseNormGMPRotation(x, rows);
zeroRows = abs(x(rows)) == 0;
assert(all(rotation(zeroRows) == 1));
assert(norm(rotation(~zeroRows) - ...
    conj(x(rows(~zeroRows)))./abs(x(rows(~zeroRows))), Inf) <= 1e-13);

%% Complex GMP and independent PN-IQ expose the same selected regressors
manager = GMP_createRegressorManager(x, y, cfg.gmp);
assert(numel(manager.regPopulation) == 673);
support = (1:12).';
U = buildGMPRegressorRows(x, rows, manager, support);
complexFit = fitComplexGMPGrid(U, y(rows), support, 1e-4);
complexPrediction = U*complexFit.coefficients;
assert(all(isfinite(complexPrediction)));

[rawFeatures, details] = buildPhaseNormalizedIQRegressors( ...
    x, rows, manager, support);
[features, reduction] = removeStructurallyZeroQFeatures( ...
    rawFeatures, details.featureMetadata, 1e-12);
normalizedTarget = rotation.*y(rows);
featureCount = min(8, size(features, 2));
selection = selectSharedIQFeatures( ...
    features, normalizedTarget, featureCount, cfg.gmp.dompOptions);
selectedFeatures = features(:, selection.supportFeatures);
pnFit = fitIndependentIQGMP( ...
    selectedFeatures, normalizedTarget, 1e-4, reduction, ...
    "Scientific smoke PN-IQ");
normalizedPrediction = selectedFeatures*pnFit.coefficientsI + ...
    1j*(selectedFeatures*pnFit.coefficientsQ);
pnPrediction = conj(rotation).*normalizedPrediction;
assert(all(isfinite(pnPrediction)));
assert(isfinite(nmseComplexDb(y(rows), pnPrediction)));
assert(details.maxCanonicalIError <= 1e-12);
assert(details.maxCanonicalQError <= 1e-12);

%% PNNN uses one visible 84-feature, two-output phase-normalized dataset
[neuralFeatures, neuralTargets, neuralRotation] = ...
    buildPhaseNormDataset(x, y, cfg.pnnn.M, ...
    cfg.pnnn.orders, cfg.pnnn.featMode);
assert(isequal(size(neuralFeatures.'), [n 84]));
assert(isequal(size(neuralTargets.'), [n 2]));
assert(numel(neuralRotation) == n);

fprintf('SCIENTIFIC PIPELINE SMOKE TEST: PASS\n');
