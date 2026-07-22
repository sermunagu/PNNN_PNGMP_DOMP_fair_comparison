% Verify unit-column, unit-peak maxima for linear families and native PNNN.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn'));

normalizedComplex = [3 + 4i; -8 + 2i; 1 - 7i];
outputPeak = 2;
comparisonComplex = normalizedComplex / outputPeak;
complexMaximum = max([abs(real(comparisonComplex)); ...
    abs(imag(comparisonComplex))]);
assert(complexMaximum == 4);
assert(abs(complexMaximum - max(abs(comparisonComplex))) > 0.01);

normalizedI = [-2; 6; 1];
normalizedQ = [4; -9; 3];
comparisonI = normalizedI / outputPeak;
comparisonQ = normalizedQ / outputPeak;
pnMaximum = max(abs([comparisonI; comparisonQ]));
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
    'fit_pniq_gmp.m'));
fixedSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'run_fixed_ridge_sweep.m'));
assert(contains(complexSource, 'normalizedCoefficients / outputPeak'));
assert(contains(complexSource, 'abs(real(activeCoefficients))'));
assert(contains(complexSource, 'abs(imag(activeCoefficients))'));
assert(contains(pnSource, 'normalizedI / outputPeak'));
assert(contains(pnSource, 'normalizedQ / outputPeak'));
assert(contains(fixedSource, 'normalizedCoefficients / outputPeak'));

fprintf('MAX ABS REAL PARAMETER TEST: PASS\n');
