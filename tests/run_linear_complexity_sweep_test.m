% Test nested Complex GMP and genuine PN-DOMP paths on a small fixture.
% The fixture verifies exact prefixes and validation-only lambda selection.
% Counters prove that matrices and maximum paths are shared across targets.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
addpath(fullfile(projectRoot, 'toolbox', 'complexity'));
addpath(fullfile(projectRoot, 'toolbox', 'domp'));
addpath(fullfile(projectRoot, 'toolbox', 'metrics'));
addpath(fullfile(projectRoot, 'toolbox', 'pn_gmp_comparison'));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

rng(913, 'twister');
n = 640;
x = complex(randn(n, 1), randn(n, 1));
y = 0.7*x + 0.12*x.*abs(x).^2 + ...
    0.03*circshift(x, 1).*abs(circshift(x, 2)).^2;
cfg = getFairDOMPComparisonConfig(projectRoot);
cfg.sweep.parameterGrid = [4 6 8];
cfg.sweep.candidateBlockSize = 128;
cfg.gmp.blockSize = 128;
split.internalTrainIndices = (1:320).';
split.internalValidationIndices = (321:400).';
split.identificationIndices = (1:480).';
split.fullSignalIndices = (1:n).';

sweep = run_linear_sweep(x, y, split, cfg);
assert(sweep.coefficientRangeDefinition == ...
    cfg.sweep.coefficientRangeDefinition);
assert(isequal(sweep.complexTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(isequal(sweep.pnTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(all(ismember(sweep.complexTable.SelectedLambda, cfg.lambdaGrid)));
assert(all(ismember(sweep.pnTable.SelectedLambda, cfg.lambdaGrid)));
assert(numel(sweep.paths.complex) == max(cfg.sweep.parameterGrid)/2);
assert(numel(sweep.paths.pn) == max(cfg.sweep.parameterGrid)/2);
assert(numel(unique(sweep.paths.complex)) == numel(sweep.paths.complex));
assert(numel(unique(sweep.paths.pn)) == numel(sweep.paths.pn));
assert(all(isfinite(sweep.predictions.complexFull), 'all'));
assert(all(isfinite(sweep.predictions.pnFull), 'all'));
assert(all(isfinite(sweep.complexTable.MaxAbsRealParameter)));
assert(all(isfinite(sweep.pnTable.MaxAbsRealParameter)));
assert(all(sweep.complexTable.MaxAbsRealParameter >= 0));
assert(all(sweep.pnTable.MaxAbsRealParameter >= 0));
assert(all(ismember({'SourceRegressorIndex','IsQ'}, ...
    sweep.pnPathMap.Properties.VariableNames)));
assert(height(sweep.pnPathMap) == numel(sweep.paths.pn));

%% Equivalent coefficients are invariant to global input/output RMS scaling
inputScale = 1.7;
outputScale = 0.6;
scaledSweep = run_linear_sweep(inputScale*x, outputScale*y, split, cfg);
assert(isequal(scaledSweep.paths, sweep.paths));
assert(isequal(scaledSweep.complexTable.SelectedLambda, ...
    sweep.complexTable.SelectedLambda));
assert(isequal(scaledSweep.pnTable.SelectedLambda, ...
    sweep.pnTable.SelectedLambda));
invariantColumns = {'ActualRealParameters','FLOPsPerSample'};
assert(isequal(scaledSweep.complexTable(:, invariantColumns), ...
    sweep.complexTable(:, invariantColumns)));
assert(isequal(scaledSweep.pnTable(:, invariantColumns), ...
    sweep.pnTable(:, invariantColumns)));
nmseColumns = {'InternalValidationNMSEdB','IdentificationNMSEdB', ...
    'FullSignalNMSEdB'};
complexNMSE = sweep.complexTable{:, nmseColumns};
scaledComplexNMSE = scaledSweep.complexTable{:, nmseColumns};
assert(all(abs(scaledComplexNMSE - complexNMSE) < 1e-8 | ...
    (scaledComplexNMSE < -250 & complexNMSE < -250), 'all'));
assert(all(abs(scaledSweep.pnTable{:, nmseColumns} - ...
    sweep.pnTable{:, nmseColumns}) < 1e-8, 'all'));
assert(all(abs(scaledSweep.complexTable.MaxAbsRealParameter - ...
    sweep.complexTable.MaxAbsRealParameter) < 1e-8));
assert(all(abs(scaledSweep.pnTable.MaxAbsRealParameter - ...
    sweep.pnTable.MaxAbsRealParameter) < 1e-8));
assert(all(abs(scaledSweep.predictions.complexFull - ...
    outputScale*sweep.predictions.complexFull) < 1e-8, 'all'));
assert(all(abs(scaledSweep.predictions.pnFull - ...
    outputScale*sweep.predictions.pnFull) < 1e-8, 'all'));

fprintf('LINEAR COMPLEXITY SWEEP TEST: PASS\n');
