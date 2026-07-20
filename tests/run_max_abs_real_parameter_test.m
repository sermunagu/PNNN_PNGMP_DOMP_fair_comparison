% Verify unit-RMS equivalent maxima for linear families and native PNNN.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn', 'pruning'));

complexCoefficients = [3 + 4i; -8 + 2i; 1 - 7i];
inputRMS = 2;
outputRMS = 8;
degrees = [1; 2; 3];
coefficientScales = inputRMS.^degrees/outputRMS;
equivalentComplex = complexCoefficients .* coefficientScales;
complexMaximum = max([abs(real(equivalentComplex)); ...
    abs(imag(equivalentComplex))]);
assert(complexMaximum == 7);
assert(abs(complexMaximum - max(abs(equivalentComplex))) > 0.01);

coefficientsI = [-2; 6; 1];
coefficientsQ = [4; -9; 3];
equivalentI = coefficientsI .* coefficientScales;
equivalentQ = coefficientsQ .* coefficientScales;
pnMaximum = max(abs([equivalentI; equivalentQ]));
assert(pnMaximum == 4.5);

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

complexSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_complex_gmp_domp.m'));
pnSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'fit_independent_pniq_domp.m'));
assert(contains(complexSource, 'inputRMS^degree/outputRMS'));
assert(contains(complexSource, 'abs(real(equivalentCoefficients))'));
assert(contains(complexSource, 'abs(imag(equivalentCoefficients))'));
assert(contains(pnSource, 'selectedMetadata.SourceRegressorIndex'));
assert(contains(pnSource, 'equivalentCoefficientsI'));
assert(contains(pnSource, 'equivalentCoefficientsQ'));

fprintf('MAX ABS REAL PARAMETER TEST: PASS\n');
