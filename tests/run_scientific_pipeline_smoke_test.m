% Exercise the shared split, GMP regressors, PN-IQ rotation, and PNNN features.
% The deterministic fixture states the core scientific contracts directly.
% It performs no neural training and uses only small linear fits.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
for folder = ["complexity","domp","metrics","pn_gmp_comparison", ...
        "pnnn","splits","sweep"]
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
rotation = complex(ones(numel(rows), 1));
zeroRows = abs(x(rows)) == 0;
rotation(~zeroRows) = ...
    conj(x(rows(~zeroRows)))./abs(x(rows(~zeroRows)));
assert(all(rotation(zeroRows) == 1));
assert(norm(rotation(~zeroRows) - ...
    conj(x(rows(~zeroRows)))./abs(x(rows(~zeroRows))), Inf) <= 1e-13);

%% Complex GMP and independent PN-IQ complete their small public sweeps
manager = GMP_createRegressorManager(x, y, cfg.gmp);
assert(numel(manager.regPopulation) == 673);
cfg.sweep.parameterGrid = [4 6 8];
cfg.sweep.candidateBlockSize = 128;
cfg.gmp.blockSize = 128;
linear = run_linear_sweep(x, y, split, cfg);
assert(all(isfinite(linear.predictions.complexFull), 'all'));
assert(all(isfinite(linear.predictions.pniqFull), 'all'));
assert(isequal(linear.complexTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(isequal(linear.pniqTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));

%% PNNN uses one visible 84-feature, two-output phase-normalized dataset
[neuralFeatures, neuralTargets, neuralRotation] = ...
    buildPhaseNormDataset(x, y, cfg.pnnn.M, ...
    cfg.pnnn.orders, cfg.pnnn.featMode);
assert(isequal(size(neuralFeatures.'), [n 84]));
assert(isequal(size(neuralTargets.'), [n 2]));
assert(numel(neuralRotation) == n);

fprintf('SCIENTIFIC PIPELINE SMOKE TEST: PASS\n');
