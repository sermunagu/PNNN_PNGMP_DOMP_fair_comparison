% Verify the common centered Welch calculation used by selected-point plots.
% Main and fixed-ridge predictions share one target reference and frequency grid.
% The fixture performs no fitting, DOMP selection, or PNNN training.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

rng(817, 'twister');
n = 1024;
sampleRateHz = 614.4e6;
time = (0:n-1).'/sampleRateHz;
target = exp(1j*2*pi*30e6*time) + ...
    0.35*exp(-1j*2*pi*75e6*time);
mainPredictions = [ ...
    target + 0.02*(randn(n, 1) + 1j*randn(n, 1)), ...
    0.99*target + 0.01*(randn(n, 1) + 1j*randn(n, 1)), ...
    0.94*target + 0.04*(randn(n, 1) + 1j*randn(n, 1))];
fixedPredictions = mainPredictions(:, [1 1 1 2 2 2]) + ...
    (1:6)*1e-3.*(randn(n, 1) + 1j*randn(n, 1));

spectrum = computeSelectedPointSpectra( ...
    target, mainPredictions, sampleRateHz, fixedPredictions);
expectedPSD = 10*log10(max(spectrum.psdLinear, ...
    spectrum.numericalFloor)/spectrum.referencePSD);

assert(size(spectrum.outputPSDdB, 2) == 4);
assert(size(spectrum.fixedOutputPSDdB, 2) == 6);
assert(size(spectrum.errorPSDdB, 2) == 3);
assert(size(spectrum.fixedErrorPSDdB, 2) == 6);
assert(isequal(spectrum.errors, target - mainPredictions));
assert(isequal(spectrum.fixedErrors, target - fixedPredictions));
assert(norm(spectrum.outputPSDdB - expectedPSD(:, 1:4), Inf) <= 1e-12);
assert(norm(spectrum.fixedOutputPSDdB - ...
    expectedPSD(:, 5:10), Inf) <= 1e-12);
assert(norm(spectrum.errorPSDdB - expectedPSD(:, 11:13), Inf) <= 1e-12);
assert(norm(spectrum.fixedErrorPSDdB - ...
    expectedPSD(:, 14:19), Inf) <= 1e-12);
assert(abs(max(spectrum.outputPSDdB(:, 1))) <= 1e-12);
assert(spectrum.config.sampleRateHz == sampleRateHz);
assert(spectrum.config.windowLength == n);
assert(spectrum.config.overlapLength == n/2);
assert(spectrum.config.nfft == 4096);
assert(spectrum.config.normalizationReference == ...
    "Maximum target full-signal PSD");

fprintf('SELECTED POINT SPECTRUM TEST: PASS\n');
