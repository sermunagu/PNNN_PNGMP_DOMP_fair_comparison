% Validate IEEE colors, NMSE limits, layouts, data references, and formats.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));
addpath(fullfile(projectRoot, 'third_party', 'matlab2tikz', 'src'));
style = getIEEEPaperStyle();
assert(isequal(style.targetGray, [117 120 123]/255));
assert(isequal(style.gmpBlue, [0 98 155]/255));
assert(isequal(style.pnOrange, [232 119 34]/255));
assert(isequal(style.pnnnGreen, [0 132 61]/255));
assert(isequal(style.selectedRed, [186 12 47]/255));

outputDirectory = getenv('PAPER_FIGURE_TEST_OUTPUT');
removeOutput = isempty(outputDirectory);
if removeOutput
    outputDirectory = tempname;
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
cleanup = onCleanup(@() removeFixture(outputDirectory, removeOutput));
budgets = [100; 200; 300; 400];
Model = [repmat("Complex GMP DOMP sweep", 4, 1); ...
    repmat("Independent PN-IQ PN-DOMP sweep", 4, 1); ...
    repmat("Sparse PNNN N12", 4, 1)];
ActualRealParameters = repmat(budgets, 3, 1);
FullSignalNMSEdB = [-34; -35; -36; -37; ...
    -33.5; -35.2; -36.5; -37.4; -34.1; -36; -38; -38.1];
FLOPsPerSample = [400; 500; 600; 700; ...
    350; 450; 550; 650; 300; 400; 500; 600];
MaxAbsRealParameter = 0.02 + ActualRealParameters/1000;
results = table(Model, ActualRealParameters, FullSignalNMSEdB, ...
    FLOPsPerSample, MaxAbsRealParameter);

lambdas = [1e-3 1e-4 1e-5];
variantBudgets = repelem(budgets, numel(lambdas));
Model = [repmat("Complex GMP-DOMP", numel(variantBudgets), 1); ...
    repmat("PN-IQ PN-DOMP", numel(variantBudgets), 1)];
ActualRealParameters = [variantBudgets; variantBudgets];
FixedLambda = repmat(repmat(lambdas(:), numel(budgets), 1), 2, 1);
FullSignalNMSEdB = -34 - ActualRealParameters/120 + ...
    0.05*repmat((1:numel(lambdas)).', 2*numel(budgets), 1);
MaxAbsRealParameter = 0.01 + ActualRealParameters/800;
fixed = table(Model, ActualRealParameters, FixedLambda, ...
    FullSignalNMSEdB, MaxAbsRealParameter);

selection = selectOperatingPoint(results, struct( ...
    'stabilizationWindowParameters', 100, ...
    'stabilizationToleranceDb', 1.00, ...
    'sensitivityWindowsParameters', [100 200], ...
    'sensitivityTolerancesDb', [0.50 1.00]));
options = struct('metricVariable', 'FullSignalNMSEdB', ...
    'metricLabel', 'Full-signal NMSE (dB)', 'includeFixed', true, ...
    'fixedLambdas', lambdas, 'isNMSE', true, ...
    'annotateSelected', false, 'selection', selection);
sweepFiles = plotSweepPaperFigure(results, fixed, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(outputDirectory, 'nmse_contract'), options);
assertFourFormats(sweepFiles);

figureHandle = openfig(char(sweepFiles.fig), 'invisible');
figureCleanup = onCleanup(@() closeFigureIfValid(figureHandle));
axesHandle = findPaperAxes(figureHandle, 'Full-signal NMSE (dB)');
limits = ylim(axesHandle);
assert(limits(2) == -30);
assert(limits(1) <= floor(min(fixed.FullSignalNMSEdB)/5)*5);
assertLineColor(axesHandle, 'Complex GMP-DOMP', style.gmpBlue);
assertLineColor(axesHandle, 'PN-IQ PN-DOMP', style.pnOrange);
assertLineColor(axesHandle, 'Sparse PNNN N12', style.pnnnGreen);
selectionMarkers = findall(axesHandle, 'Type', 'line', ...
    'DisplayName', 'All three models used in selection');
assert(isscalar(selectionMarkers));
assert(numel(selectionMarkers.XData) == 3);
assert(max(abs(selectionMarkers.Color - style.selectedRed)) < 1e-12);
pnHighlight = findall(axesHandle, 'Type', 'line', ...
    'DisplayName', 'PN-IQ at selected common budget');
assert(isscalar(pnHighlight));
assert(max(abs(pnHighlight.Color - style.selectedRed)) < 1e-12);
assert(strcmpi(pnHighlight.HandleVisibility, 'off'));
legendHandle = findobj(figureHandle, 'Type', 'legend');
assert(isscalar(legendHandle));
assert(strcmpi(legendHandle.Orientation, 'horizontal'));
clear figureCleanup;

frequencyMHz = linspace(-100, 100, 128).';
targetPSD = -abs(frequencyMHz)/8;
spectrum = struct( ...
    'frequencyMHz', frequencyMHz, ...
    'outputPSDdB', [targetPSD, targetPSD-1, targetPSD-2, targetPSD-3], ...
    'fixedOutputPSDdB', targetPSD - (1:6), ...
    'errorPSDdB', -45 - abs(frequencyMHz)/12 - (1:3), ...
    'fixedErrorPSDdB', -47 - abs(frequencyMHz)/13 - (1:6));
spectrumFiles = exportSelectedSpectrumFigures( ...
    spectrum, outputDirectory, lambdas, struct());
fields = fieldnames(spectrumFiles);
for index = 1:numel(fields)
    assertFourFormats(spectrumFiles.(fields{index}));
end

figureHandle = openfig(char(spectrumFiles.error.fig), 'invisible');
figureCleanup = onCleanup(@() closeFigureIfValid(figureHandle));
axesHandle = findPaperAxes(figureHandle, ...
    'PSD relative to target peak (dB)');
dataLines = findobj(axesHandle, 'Type', 'line');
assert(numel(dataLines) == 4);
assertLineColor(axesHandle, 'Target full signal', style.targetGray);
legendHandle = findobj(figureHandle, 'Type', 'legend');
assert(strcmpi(legendHandle.Orientation, 'horizontal'));
clear figureCleanup;

figureHandle = openfig(char(spectrumFiles.ridgeError.fig), 'invisible');
figureCleanup = onCleanup(@() closeFigureIfValid(figureHandle));
axesHandles = findobj(figureHandle, 'Type', 'axes');
assert(numel(axesHandles) == 2);
positions = vertcat(axesHandles.Position);
assert(abs(diff(positions(:, 2))) < 0.1);
assert(abs(diff(positions(:, 1))) > 0.2);
for index = 1:2
    assertLineColor(axesHandles(index), 'Target full signal', ...
        style.targetGray);
end
legendHandle = findobj(figureHandle, 'Type', 'legend');
assert(isscalar(legendHandle));
assert(strcmpi(legendHandle.Orientation, 'horizontal'));
clear figureCleanup;

clear cleanup;
fprintf('PAPER FIGURE CONTRACT TEST: PASS\n');

function axesHandle = findPaperAxes(figureHandle, yLabel)
axesHandles = findobj(figureHandle, 'Type', 'axes');
matches = false(size(axesHandles));
for index = 1:numel(axesHandles)
    matches(index) = string(axesHandles(index).YLabel.String) == yLabel;
end
axesHandle = axesHandles(matches);
assert(isscalar(axesHandle));
end


function closeFigureIfValid(figureHandle)
if isgraphics(figureHandle, 'figure')
    close(figureHandle);
end
end

function assertLineColor(axesHandle, displayName, expectedColor)
lineHandle = findobj(axesHandle, 'Type', 'line', ...
    'DisplayName', displayName);
assert(isscalar(lineHandle));
assert(max(abs(lineHandle.Color - expectedColor)) < 1e-12);
end

function assertFourFormats(files)
assert(isfile(files.fig));
assert(isfile(files.png));
assert(isfile(files.tikz));
assert(isfile(files.pdf));
end

function removeFixture(directory, removeOutput)
close all force;
if removeOutput && isfolder(directory)
    rmdir(directory, 's');
end
end
