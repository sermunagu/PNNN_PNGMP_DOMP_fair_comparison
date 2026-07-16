% Test nested Complex GMP and genuine PN-DOMP paths on a small fixture.
% The fixture verifies exact prefixes and validation-only lambda selection.
% A controlled support separately checks the fixed historical DOMP-100 contract.

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
cfg.sweep.includeHistoricalPNIQReference = false;
cfg.sweep.candidateBlockSize = 128;
cfg.gmp.blockSize = 128;
split.internalTrainIndices = (1:320).';
split.internalValidationIndices = (321:400).';
split.identificationIndices = (1:480).';
split.fullSignalIndices = (1:n).';

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
assert(isempty(sweep.historicalTable));

manager = GMP_createRegressorManager(x, y, cfg.gmp);
population = (1:numel(manager.regPopulation)).';
[~, structure] = analyzeRegressorStructure( ...
    manager.regPopulation, population);
descriptors = structure.descriptors;
realIndices = find([descriptors.QColumnStructurallyZero]);
complexIndices = find(~[descriptors.QColumnStructurallyZero]);
realIndices = realIndices(:);
complexIndices = complexIndices(:);
historicalSupport = [realIndices(1:28); complexIndices(1:72)];
[historicalRaw, historicalDetails] = buildPhaseNormalizedIQRegressors( ...
    x, split.internalTrainIndices, manager, historicalSupport);
[historicalFeatures, ~] = removeStructurallyZeroQFeatures( ...
    historicalRaw, historicalDetails.featureMetadata, 1e-12);
assert(cfg.gmp.baseSupportSize == 100);
assert(size(historicalRaw, 2) == 200);
assert(size(historicalFeatures, 2) == 172);

remaining = setdiff(population, historicalSupport, 'stable');
maximumPath = [historicalSupport; remaining];
grids = {[300 344 380], [280 344 400]};
historicalSupports = cell(size(grids));
for index = 1:numel(grids)
    gridConfig = cfg;
    gridConfig.sweep.parameterGrid = grids{index};
    maximumTerms = max(max(gridConfig.sweep.parameterGrid)/2, ...
        gridConfig.gmp.baseSupportSize);
    assert(numel(maximumPath) >= maximumTerms);
    historicalSupports{index} = ...
        maximumPath(1:gridConfig.gmp.baseSupportSize);
end
[identificationRaw, identificationDetails] = ...
    buildPhaseNormalizedIQRegressors(x, split.identificationIndices, ...
    manager, historicalSupports{1});
[identificationFeatures, ~] = removeStructurallyZeroQFeatures( ...
    identificationRaw, identificationDetails.featureMetadata, 1e-12);
metadata = struct( ...
    'HistoricalGMPComplexSupportSize', numel(historicalSupports{1}), ...
    'HistoricalEffectivePNFeatures', size(identificationFeatures, 2), ...
    'HistoricalNumRealParameters', 2*size(identificationFeatures, 2));
assert(metadata.HistoricalGMPComplexSupportSize == 100);
assert(metadata.HistoricalEffectivePNFeatures == 172);
assert(metadata.HistoricalNumRealParameters == 344);
assert(isequal(historicalSupports{1}, historicalSupport));
assert(isequal(historicalSupports{1}, historicalSupports{2}));
parameterMatchedSupport = maximumPath(1:344/2);
assert(~isequal(historicalSupports{1}, parameterMatchedSupport));

fprintf('LINEAR COMPLEXITY SWEEP TEST: PASS\n');
