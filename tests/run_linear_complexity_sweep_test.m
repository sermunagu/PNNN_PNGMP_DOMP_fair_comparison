% Test one identification DOMP path and lambda-zero prefix fits per family.
% The fixture deliberately omits internal split fields from the linear API.

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
split.identificationIndices = (1:480).';
split.fullSignalIndices = (1:n).';

sweep = run_linear_sweep(x, y, split, cfg);
assert(sweep.coefficientRangeDefinition == ...
    cfg.sweep.coefficientRangeDefinition);
assert(isequal(sweep.complexTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(isequal(sweep.pniqTable.ActualRealParameters.', ...
    cfg.sweep.parameterGrid));
assert(all(sweep.complexTable.SelectedLambda == 0));
assert(all(sweep.pniqTable.SelectedLambda == 0));
assert(all(isnan(sweep.complexTable.InternalValidationNMSEdB)));
assert(all(isnan(sweep.pniqTable.InternalValidationNMSEdB)));
assert(numel(sweep.paths.complex) == max(cfg.sweep.parameterGrid)/2);
assert(numel(sweep.paths.pniq) == max(cfg.sweep.parameterGrid)/2);
assert(numel(unique(sweep.paths.complex)) == numel(sweep.paths.complex));
assert(numel(unique(sweep.paths.pniq)) == numel(sweep.paths.pniq));
assert(all(isfinite(sweep.predictions.complexFull), 'all'));
assert(all(isfinite(sweep.predictions.pniqFull), 'all'));
assert(all(isfinite(sweep.complexTable.MaxAbsRealParameter)));
assert(all(isfinite(sweep.pniqTable.MaxAbsRealParameter)));
assert(all(sweep.complexTable.MaxAbsRealParameter >= 0));
assert(all(sweep.pniqTable.MaxAbsRealParameter >= 0));
assert(all(ismember({'SourceRegressorIndex','IsQ'}, ...
    sweep.pniqPathMap.Properties.VariableNames)));
assert(height(sweep.pniqPathMap) == numel(sweep.paths.pniq));
complexPrefixes = cell(numel(cfg.sweep.parameterGrid), 1);
pniqPrefixes = cell(numel(cfg.sweep.parameterGrid), 1);
for targetIndex = 1:numel(cfg.sweep.parameterGrid)
    count = cfg.sweep.parameterGrid(targetIndex)/2;
    complexPrefixes{targetIndex} = sweep.paths.complex(1:count);
    pniqPrefixes{targetIndex} = sweep.paths.pniq(1:count);
end
for targetIndex = 2:numel(cfg.sweep.parameterGrid)
    previousCount = numel(complexPrefixes{targetIndex - 1});
    assert(isequal(complexPrefixes{targetIndex - 1}, ...
        complexPrefixes{targetIndex}(1:previousCount)));
    assert(isequal(pniqPrefixes{targetIndex - 1}, ...
        pniqPrefixes{targetIndex}(1:previousCount)));
end

%% Source contract: one DOMP call, no internal split or lambda-grid access
complexSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_complex_gmp_domp.m'));
pniqSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_pniq_gmp.m'));
for source = string({complexSource, pniqSource})
    sourceText = char(source);
    assert(isscalar(strfind(sourceText, 'selectDOMPSupport')));
    assert(~contains(source, 'internalTrainIndices'));
    assert(~contains(source, 'internalValidationIndices'));
    assert(~contains(source, 'cfg.lambdaGrid'));
    assert(contains(source, 'lsqminnorm'));
end

%% Coefficient ranges match an explicit per-column peak construction
[explicitComplex, explicitPNIQ] = explicitCoefficientRanges( ...
    x, y, split, cfg, sweep);
assert(all(abs(explicitComplex - ...
    sweep.complexTable.MaxAbsRealParameter) < 1e-9));
assert(all(abs(explicitPNIQ - ...
    sweep.pniqTable.MaxAbsRealParameter) < 1e-9));

%% Equivalent normalized coefficients are invariant to input/output scaling
inputScale = 1.7;
outputScale = 0.6;
scaledSweep = run_linear_sweep(inputScale*x, outputScale*y, split, cfg);
assert(isequal(scaledSweep.paths, sweep.paths));
assert(isequal(scaledSweep.complexTable.SelectedLambda, ...
    sweep.complexTable.SelectedLambda));
assert(isequal(scaledSweep.pniqTable.SelectedLambda, ...
    sweep.pniqTable.SelectedLambda));
invariantColumns = {'ActualRealParameters','FLOPsPerSample'};
assert(isequal(scaledSweep.complexTable(:, invariantColumns), ...
    sweep.complexTable(:, invariantColumns)));
assert(isequal(scaledSweep.pniqTable(:, invariantColumns), ...
    sweep.pniqTable(:, invariantColumns)));
nmseColumns = {'IdentificationNMSEdB','FullSignalNMSEdB'};
complexNMSE = sweep.complexTable{:, nmseColumns};
scaledComplexNMSE = scaledSweep.complexTable{:, nmseColumns};
assert(all(abs(scaledComplexNMSE - complexNMSE) < 1e-8 | ...
    (scaledComplexNMSE < -250 & complexNMSE < -250), 'all'));
assert(all(abs(scaledSweep.pniqTable{:, nmseColumns} - ...
    sweep.pniqTable{:, nmseColumns}) < 1e-8, 'all'));
assert(all(abs(scaledSweep.complexTable.MaxAbsRealParameter - ...
    sweep.complexTable.MaxAbsRealParameter) < 1e-8));
assert(all(abs(scaledSweep.pniqTable.MaxAbsRealParameter - ...
    sweep.pniqTable.MaxAbsRealParameter) < 1e-8));
assert(all(abs(scaledSweep.predictions.complexFull - ...
    outputScale*sweep.predictions.complexFull) < 1e-8, 'all'));
assert(all(abs(scaledSweep.predictions.pniqFull - ...
    outputScale*sweep.predictions.pniqFull) < 1e-8, 'all'));

fprintf('LINEAR COMPLEXITY SWEEP TEST: PASS\n');

function [complexRanges, pniqRanges] = explicitCoefficientRanges( ...
    x, y, split, cfg, sweep)
rows = split.identificationIndices(:);
outputPeak = max(abs(y(rows)));
unitPeakTarget = y / outputPeak;
manager = GMP_createRegressorManager(x, y, cfg.gmp);
targets = cfg.sweep.parameterGrid(:);
complexRanges = zeros(numel(targets), 1);
pniqRanges = zeros(numel(targets), 1);

for targetIndex = 1:numel(targets)
    count = targets(targetIndex)/2;
    support = sweep.paths.complex(1:count);
    regressors = buildGMPRegressorRows( ...
        x, rows, manager, support);
    regressors = regressors ./ max(abs(regressors), [], 1);
    coefficients = fitUnitColumns(regressors, unitPeakTarget(rows));
    complexRanges(targetIndex) = max(abs(coefficients));

    metadata = sweep.pniqPathMap(1:count, :);
    complexSupport = unique(metadata.SourceRegressorIndex, 'stable');
    complexRegressors = buildGMPRegressorRows( ...
        x, rows, manager, complexSupport);
    rotation = complex(ones(numel(rows), 1));
    nonzero = abs(x(rows)) ~= 0;
    rotation(nonzero) = conj(x(rows(nonzero))) ./ abs(x(rows(nonzero)));
    phaseNormalized = rotation .* complexRegressors;
    [~, sourceColumns] = ismember( ...
        metadata.SourceRegressorIndex, complexSupport);
    features = zeros(numel(rows), count);
    for featureIndex = 1:count
        values = phaseNormalized(:, sourceColumns(featureIndex));
        if metadata.IsQ(featureIndex)
            features(:, featureIndex) = imag(values);
        else
            features(:, featureIndex) = real(values);
        end
    end
    features = features ./ max(abs(features), [], 1);
    target = rotation .* unitPeakTarget(rows);
    coefficientsI = fitUnitColumns(features, real(target));
    coefficientsQ = fitUnitColumns(features, imag(target));
    pniqRanges(targetIndex) = max(abs([coefficientsI; coefficientsQ]));
end
end

function coefficients = fitUnitColumns(regressors, target)
tolerance = max(size(regressors))*eps(norm(regressors, 2));
coefficients = lsqminnorm(regressors, target, tolerance);
end
