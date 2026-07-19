% Test the canonical sweep table and presentation contract.
% The 49-target fixture verifies three principal families and six fixed-Ridge
% variants without loading measurements, fitting models, or opening figures.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'config'));
cfg = getFairDOMPComparisonConfig(projectRoot);
targets = (20:10:500).';
targetCount = numel(targets);

Model = [repmat("Complex GMP DOMP sweep", targetCount, 1); ...
    repmat("Independent PN-IQ PN-DOMP sweep", targetCount, 1); ...
    repmat("Sparse PNNN N12", targetCount, 1)];
TargetRealParameters = repmat(targets, 3, 1);
ActualRealParameters = TargetRealParameters;
FullSignalNMSEdB = [-25 - 0.01*targets; ...
    -25.5 - 0.01*targets; -22 - 0.008*targets];
FLOPsPerSample = [700 + targets; 500 + targets; 400 + targets];
ActiveWeights = ActualRealParameters;
ActiveBiases = zeros(3*targetCount, 1);
WeightSparsityPercent = nan(3*targetCount, 1);
SelectedLambda = nan(3*targetCount, 1);
InternalValidationNMSEdB = nan(3*targetCount, 1);
IdentificationNMSEdB = FullSignalNMSEdB + 0.1;
FineTuneEpochs = nan(3*targetCount, 1);
MaxAbsRealParameter = 0.01 + ActualRealParameters/1000;
pnnnRows = Model == "Sparse PNNN N12";
ActiveBiases(pnnnRows) = 14;
ActiveWeights(pnnnRows) = ActualRealParameters(pnnnRows) - 14;
results = table(Model, TargetRealParameters, ...
    ActualRealParameters, SelectedLambda, InternalValidationNMSEdB, ...
    IdentificationNMSEdB, FullSignalNMSEdB, FLOPsPerSample, ...
    ActiveWeights, ActiveBiases, WeightSparsityPercent, FineTuneEpochs, ...
    MaxAbsRealParameter);

fixedLambdas = cfg.fixedRidgeLambdas(:);
variantTargets = repelem(targets, numel(fixedLambdas));
Model = [repmat("Complex GMP-DOMP", numel(variantTargets), 1); ...
    repmat("PN-IQ PN-DOMP", numel(variantTargets), 1)];
TargetRealParameters = [variantTargets; variantTargets];
ActualRealParameters = TargetRealParameters;
FixedLambda = repmat(repmat(fixedLambdas, targetCount, 1), 2, 1);
FLOPsPerSample = [repelem(results.FLOPsPerSample(1:targetCount), 3); ...
    repelem(results.FLOPsPerSample(targetCount + (1:targetCount)), 3)];
IdentificationNMSEdB = -35 - 0.001*TargetRealParameters;
FullSignalNMSEdB = IdentificationNMSEdB - 0.1;
MaxAbsRealParameter = 0.02 + ActualRealParameters/900;
fixedResults = table(Model, TargetRealParameters, ...
    ActualRealParameters, FixedLambda, IdentificationNMSEdB, ...
    FullSignalNMSEdB, FLOPsPerSample, MaxAbsRealParameter);

assert(height(results) == 147);
assert(height(fixedResults) == 294);
assert(width(results) == 13);
assert(width(fixedResults) == 8);
assert(isequal(sort(unique(results.Model)), sort([ ...
    "Complex GMP DOMP sweep"; "Independent PN-IQ PN-DOMP sweep"; ...
    "Sparse PNNN N12"])));
assert(~any(contains(results.Model, "Historical")));
for target = targets.'
    rows = results.TargetRealParameters == target;
    assert(nnz(rows) == 3);
    assert(all(results.ActualRealParameters(rows) == target));
end
assert(all(results.ActiveWeights(pnnnRows) + ...
    results.ActiveBiases(pnnnRows) == results.ActualRealParameters(pnnnRows)));
assert(all(isnan(results.WeightSparsityPercent(~pnnnRows))));
for model = unique(fixedResults.Model).'
    for lambda = fixedLambdas.'
        rows = fixedResults.Model == model & ...
            fixedResults.FixedLambda == lambda;
        assert(nnz(rows) == targetCount);
    end
end
assert(all(results.MaxAbsRealParameter >= 0));
assert(all(fixedResults.MaxAbsRealParameter >= 0));

fprintf('SWEEP PRESENTATION TEST: PASS\n');
