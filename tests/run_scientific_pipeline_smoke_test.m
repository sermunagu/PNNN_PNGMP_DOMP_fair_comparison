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
% Fixed P200 support from legacy/sweep_4329deabf42b; regression belongs here,
% not inside the operation counter for every arbitrary 100-term model.
support = [1 12 8 23 34 24 10 11 2 437 434 3 45 459 408 173 9 56 ...
    33 503 413 410 479 166 7 6 5 4 27 476 441 89 457 454 523 145 ...
    417 520 672 431 187 44 567 320 456 57 298 588 101 569 255 430 ...
    167 465 564 154 328 442 174 28 445 175 205 67 100 321 452 518 ...
    411 36 478 254 500 505 190 458 211 461 438 26 169 170 517 189 ...
    593 13 440 171 172 571 597 78 65 155 501 68 498 443 429 186];
operations = countModelOperations(manager.regPopulation, support);
assert(isequal(operations.TotalRealMultiplicationsPerSample, [666;1074;632;632]));
assert(isequal(operations.TotalRealAdditionsPerSample, [413;617;387;387]));
assert(all(isfinite(operations{:,3:end}), 'all'));
assert(all(operations{:,3:end} >= 0, 'all'));
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
