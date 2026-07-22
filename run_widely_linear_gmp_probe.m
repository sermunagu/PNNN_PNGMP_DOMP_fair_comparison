function comparison = run_widely_linear_gmp_probe()
% Compare Complex GMP-DOMP, WL-GMP, and PN-IQ-GMP at 340 real parameters.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'toolbox')));
cfg = getFairDOMPComparisonConfig(root);

target = 340;
nAtoms = target/2;
data = load(cfg.measurementFile, 'x', 'y');
[x, y] = selectXYByMapping(data.x, data.y, cfg.mappingMode);
x = x(:); y = y(:);
if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end

split = buildCommonComparisonSplit(x, y, cfg);
manager = GMP_createRegressorManager(x, y, cfg.gmp);
population = (1:numel(manager.regPopulation)).';

trainRows = split.internalTrainIndices(:);
validationRows = split.internalValidationIndices(:);
identificationRows = split.identificationIndices(:);
fullRows = split.fullSignalIndices(:);

U = buildGMPRegressorRows(x, trainRows, manager, population);
trainU = [U, conj(U)];

U = buildGMPRegressorRows(x, validationRows, manager, population);
validationU = [U, conj(U)];

trainPath = selectDOMPSupport(trainU, y(trainRows), nAtoms, cfg.gmp.dompOptions.columnTolerance);
assert(numel(trainPath) == nAtoms, 'WL-GMP could not select 170 independent atoms.');

validationNMSE = zeros(numel(cfg.lambdaGrid), 1);
for k = 1:numel(cfg.lambdaGrid)
    c = ridgeFit(trainU(:, trainPath), y(trainRows), cfg.lambdaGrid(k));
    validationNMSE(k) = nmseComplexDb(y(validationRows), validationU(:, trainPath)*c);
end

[~, best] = min(validationNMSE);
lambda = cfg.lambdaGrid(best);
clear trainU validationU U

U = buildGMPRegressorRows(x, identificationRows, manager, population);
identificationU = [U, conj(U)];
path = selectDOMPSupport(identificationU, y(identificationRows), nAtoms, cfg.gmp.dompOptions.columnTolerance);
assert(numel(path) == nAtoms, 'WL-GMP could not select 170 independent atoms.');

c = ridgeFit(identificationU(:, path), y(identificationRows), lambda);
identificationNMSE = nmseComplexDb(y(identificationRows), identificationU(:, path)*c);
clear identificationU U

prediction = complex(zeros(numel(fullRows), 1));
for first = 1:cfg.gmp.blockSize:numel(fullRows)
    local = first:min(first + cfg.gmp.blockSize - 1, numel(fullRows));
    U = buildGMPRegressorRows(x, fullRows(local), manager, population);
    U = [U, conj(U)]; %#ok<AGROW>
    prediction(local) = U(:, path)*c;
end

wlNMSE = nmseComplexDb(y(fullRows), prediction);

referenceFile = fullfile(cfg.sweep.resultsRoot, 'sweep_d113e389ab78', 'complexity_sweep.csv');
reference = readtable(referenceFile);
rows = reference.TargetRealParameters == target;
referenceRows = reference(rows, :);
gmp = referenceRows(1, :);
pniq = referenceRows(2, :);

Model = [cfg.names.complexGMPDOMP; "Widely linear GMP"; cfg.names.pniqGMP];
FullSignalNMSEdB = [gmp.FullSignalNMSEdB; wlNMSE; pniq.FullSignalNMSEdB];
ImprovementVsGMPdB = gmp.FullSignalNMSEdB - FullSignalNMSEdB;
comparison = table(Model, FullSignalNMSEdB, ImprovementVsGMPdB);
disp(comparison);

fprintf('WL-GMP lambda: %.3g\n', lambda);
fprintf('WL-GMP identification NMSE: %.6f dB\n', identificationNMSE);
fprintf('WL-GMP support: %d direct + %d conjugate atoms\n', nnz(path <= numel(population)), nnz(path > numel(population)));

outputFile = fullfile(fileparts(referenceFile), 'widely_linear_gmp_probe_0340.csv');
writetable(comparison, outputFile);
fprintf('Written: %s\n', outputFile);
end


function c = ridgeFit(U, y, lambda)
norms = sqrt(sum(abs(U).^2, 1)).';
U = U ./ norms.';
if lambda == 0
    tolerance = max(size(U))*eps(norm(U, 2));
    c = lsqminnorm(U, y, tolerance);
else
    c = (U'*U + lambda*eye(size(U, 2))) \ (U'*y);
end
c = c ./ norms;
end
