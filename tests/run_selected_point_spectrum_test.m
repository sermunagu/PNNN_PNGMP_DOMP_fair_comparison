% Test shared spectral calculation and the two selected-point figures.
% A small complex fixture verifies normalization, errors, curves, and overwrite.
% No measurement, DOMP, fitting, pruning, or PNNN execution is performed.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

rng(817, 'twister');
n = 1024;
sampleRateHz = 614.4e6;
time = (0:n-1).'/sampleRateHz;
target = exp(1j*2*pi*30e6*time) + ...
    0.35*exp(-1j*2*pi*75e6*time);
predictions = [target + 0.02*(randn(n, 1) + 1j*randn(n, 1)), ...
    0.99*target + 0.01*(randn(n, 1) + 1j*randn(n, 1)), ...
    0.94*target + 0.04*(randn(n, 1) + 1j*randn(n, 1))];
fixedPredictions = predictions(:, [1 1 1 2 2 2]) + ...
    (1:6)*1e-3.*(randn(n, 1) + 1j*randn(n, 1));
spectrum = computeSelectedPointSpectra( ...
    target, predictions, sampleRateHz, fixedPredictions);

assert(size(spectrum.outputPSDdB, 1) == numel(spectrum.frequencyMHz));
assert(size(spectrum.errorPSDdB, 1) == numel(spectrum.frequencyMHz));
assert(size(spectrum.outputPSDdB, 2) == 4);
assert(size(spectrum.errorPSDdB, 2) == 3);
assert(size(spectrum.fixedOutputPSDdB, 2) == 6);
assert(size(spectrum.fixedErrorPSDdB, 2) == 6);
assert(abs(max(spectrum.outputPSDdB(:, 1))) <= 1e-12);
assert(isequal(spectrum.errors, target - predictions));
assert(isequal(spectrum.fixedErrors, target - fixedPredictions));
expectedPSD = 10*log10(max(spectrum.psdLinear, ...
    spectrum.numericalFloor)/spectrum.referencePSD);
assert(norm(spectrum.outputPSDdB - expectedPSD(:, 1:4), Inf) <= 1e-12);
assert(norm(spectrum.fixedOutputPSDdB - expectedPSD(:, 5:10), Inf) <= 1e-12);
assert(norm(spectrum.errorPSDdB - expectedPSD(:, 11:13), Inf) <= 1e-12);
assert(norm(spectrum.fixedErrorPSDdB - expectedPSD(:, 14:19), Inf) <= 1e-12);

outputRoot = tempname;
mkdir(outputRoot);
outputCleanup = onCleanup(@() rmdir(outputRoot, 's'));
checkpointFiles = ["linear_sweep.mat", "sweep_dense_source.mat", ...
    "pnnn_target_0340.mat"];
checkpointSentinel = uint8(0:31);
for filename = checkpointFiles
    save(fullfile(outputRoot, filename), 'checkpointSentinel');
end
bytesBefore = arrayfun(@(name) readBinaryFile(fullfile(outputRoot, name)), ...
    checkpointFiles, 'UniformOutput', false);

Model = ["Complex GMP-DOMP"; "PN-IQ PN-DOMP"; "Sparse PNNN N12"];
comparisonTable = table(Model);
fixedModel = repelem(Model(1:2), 3);
FixedLambda = repmat([1e-3; 1e-4; 1e-5], 2, 1);
fixedLambdaComparisonTable = table(fixedModel, FixedLambda, ...
    'VariableNames', {'Model','FixedLambda'});
fixedLambdaFullSignalPredictions = struct( ...
    'complexGMP', struct('lambda1e3', fixedPredictions(:, 1), ...
        'lambda1e4', fixedPredictions(:, 2), ...
        'lambda1e5', fixedPredictions(:, 3)), ...
    'pnIQ', struct('lambda1e3', fixedPredictions(:, 4), ...
        'lambda1e4', fixedPredictions(:, 5), ...
        'lambda1e5', fixedPredictions(:, 6)));
results = struct('selectedParameters', 340, ...
    'resultDirectory', string(outputRoot), ...
    'targetFullSignal', target, 'fullSignalIndices', (1:n).', ...
    'sampleRateHz', sampleRateHz, 'sampleRateSource', "Synthetic fixture", ...
    'fullSignalPredictions', struct('complexGMP', predictions(:, 1), ...
        'pnIQ', predictions(:, 2), 'sparsePNNNN12', predictions(:, 3)), ...
    'fixedLambdaFullSignalPredictions', ...
        fixedLambdaFullSignalPredictions, ...
    'fixedLambdaComparisonTable', fixedLambdaComparisonTable, ...
    'comparisonTable', comparisonTable);
first = writeSelectedPointSpectra(results);
second = writeSelectedPointSpectra(results);

assert(first.selectedPointDirectory == second.selectedPointDirectory);
assert(isfile(second.outputSpectrumFigure));
assert(isfile(second.errorSpectrumFigure));
assert(isfile(second.ridgeOutputSpectrumFigure));
assert(isfile(second.ridgeErrorSpectrumFigure));
assert(second.spectrumConfig.outputCurveCount == 4);
assert(second.spectrumConfig.errorCurveCount == 3);
assert(isequal(second.spectrumConfig.ridgeOutputPanelCurveCounts, [5 5]));
assert(isequal(second.spectrumConfig.ridgeErrorPanelCurveCounts, [4 4]));
assert(second.spectrumConfig.sampleRateHz == sampleRateHz);
assert(second.spectrumConfig.normalizationReference == ...
    "Maximum target full-signal PSD");
assert(height(second.comparisonTable) == 3);
assert(~any(contains(second.comparisonTable.Model, ...
    ["Ridge", "H4", "dense"], 'IgnoreCase', true), 'all'));
figures = dir(fullfile(second.selectedPointDirectory, '*.png'));
assert(numel(figures) == 4);
assert(~isempty(imfinfo(second.outputSpectrumFigure)));
assert(~isempty(imfinfo(second.errorSpectrumFigure)));
assert(~isempty(imfinfo(second.ridgeOutputSpectrumFigure)));
assert(~isempty(imfinfo(second.ridgeErrorSpectrumFigure)));
bytesAfter = arrayfun(@(name) readBinaryFile(fullfile(outputRoot, name)), ...
    checkpointFiles, 'UniformOutput', false);
assert(isequal(bytesAfter, bytesBefore));

source = string(fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'writeSelectedPointSpectra.m'))) + newline + ...
    string(fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'computeSelectedPointSpectra.m')));
for forbidden = ["selectDOMPSupport", "runPNNNComparisonStudy", ...
        "runPNNNSparseSweep", "fitComplexGMPGrid"]
    assert(~contains(source, forbidden));
end
clear outputCleanup;

fprintf('SELECTED POINT SPECTRUM TEST: PASS\n');

function bytes = readBinaryFile(filename)
file = fopen(filename, 'rb');
assert(file >= 0);
cleanup = onCleanup(@() fclose(file));
bytes = fread(file, Inf, '*uint8');
clear cleanup;
end
