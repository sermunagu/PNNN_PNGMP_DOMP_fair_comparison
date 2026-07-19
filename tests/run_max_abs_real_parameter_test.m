% Verify native stored-scalar maxima for complex, PN-IQ, and masked PNNN.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'pnnn', 'pruning'));

complexCoefficients = [3 + 4i; -8 + 2i; 1 - 7i];
complexMaximum = max([abs(real(complexCoefficients)); ...
    abs(imag(complexCoefficients))]);
assert(complexMaximum == 8);
assert(abs(complexMaximum - max(abs(complexCoefficients))) > 0.1);

coefficientsI = [-2; 6; 1];
coefficientsQ = [4; -9; 3];
pnMaximum = max(abs([coefficientsI; coefficientsQ]));
assert(pnMaximum == 9);

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
assert(contains(complexSource, 'abs(real(activeCoefficients))'));
assert(contains(complexSource, 'abs(imag(activeCoefficients))'));
assert(contains(pnSource, 'coefficientsI(1:count, targetIndex)'));
assert(contains(pnSource, 'coefficientsQ(1:count, targetIndex)'));

fprintf('MAX ABS REAL PARAMETER TEST: PASS\n');
