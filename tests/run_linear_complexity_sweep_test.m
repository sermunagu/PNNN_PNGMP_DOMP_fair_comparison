% Test nested Complex GMP and genuine PN-DOMP paths on a small fixture.
% The fixture verifies exact prefixes, validation-only lambda selection,
% finite predictions, and a separate configurable historical reference.

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
cfg.gmp.maxPopulation = 2;
cfg.gmp.blockSize = 128;
split.internalTrainIndices = (1:320).';
split.internalValidationIndices = (321:400).';
split.identificationIndices = (1:480).';
split.fullSignalIndices = (1:n).';

manager = GMP_createRegressorManager(x, y, cfg.gmp);
population = (1:numel(manager.regPopulation)).';
trainU = buildGMPRegressorRows( ...
    x, split.internalTrainIndices, manager, population);
historicalSupport = selectDOMPSupport( ...
    trainU, y(split.internalTrainIndices), cfg.gmp.maxPopulation, ...
    cfg.gmp.dompOptions);
[historicalRaw, historicalDetails] = buildPhaseNormalizedIQRegressors( ...
    x, split.internalTrainIndices, manager, historicalSupport);
[historicalFeatures, ~] = removeStructurallyZeroQFeatures( ...
    historicalRaw, historicalDetails.featureMetadata, 1e-12);
cfg.sweep.historicalReferenceParameters = 2*size(historicalFeatures, 2);

sweep = runLinearComplexitySweep(x, y, split, cfg);
assert(isequal(sweep.complexTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(isequal(sweep.pnTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(all(ismember(sweep.complexTable.SelectedLambda, cfg.lambdaGrid)));
assert(all(ismember(sweep.pnTable.SelectedLambda, cfg.lambdaGrid)));
assert(sweep.metadata.retainedFeatures >= max(cfg.sweep.parameterGrid)/2);
assert(~sweep.metadata.fullSignalUsedForSelection);
assert(~sweep.metadata.fullSignalUsedForFitting);
assert(all(struct2array(sweep.metadata.dompInvocationCount) == 1));
for supports = {sweep.supports.complex, sweep.supports.pnFeatures}
    family = supports{1};
    for index = 2:numel(family)
        assert(isequal(family{index}(1:numel(family{index-1})), ...
            family{index-1}));
    end
end
assert(all(isfinite(sweep.predictions.complexIdentification), 'all'));
assert(all(isfinite(sweep.predictions.complexFull), 'all'));
assert(all(isfinite(sweep.predictions.pnIdentification), 'all'));
assert(all(isfinite(sweep.predictions.pnFull), 'all'));
assert(height(sweep.historicalTable) == 1);
assert(sweep.historicalTable.SweepRole == "Historical reference");
assert(sweep.historicalTable.ActualRealParameters == ...
    cfg.sweep.historicalReferenceParameters);
assert(all(isfinite(sweep.predictions.historicalIdentification)));
assert(all(isfinite(sweep.predictions.historicalFull)));

fprintf('LINEAR COMPLEXITY SWEEP TEST: PASS\n');
