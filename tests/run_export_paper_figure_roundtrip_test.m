% Verify that exportPaperFigure writes FIG, PNG, and layout-neutral TikZ.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));
addpath(fullfile(projectRoot, 'third_party', 'matlab2tikz', 'src'));
outputDirectory = tempname;
mkdir(outputDirectory);
cleanup = onCleanup(@() removeFixture(outputDirectory));

figureHandle = figure('Visible', 'off');
axesHandle = axes(figureHandle);
plot(axesHandle, 1:3, [1 4 2]);
files = exportPaperFigure(figureHandle, ...
    fullfile(outputDirectory, 'roundtrip'), struct());
close(figureHandle);

assert(isequal(sort(string(fieldnames(files))), sort(["fig"; "png"; "tikz"])));
assert(isfile(files.fig));
assert(isfile(files.png));
assert(isfile(files.tikz));
tikzText = fileread(files.tikz);
assert(contains(tikzText, '\figurewidth'));
assert(contains(tikzText, '\figureheight'));
assert(~contains(tikzText, '0.68\figurewidth'));
assert(~contains(tikzText, '0.7\figureheight'));
absoluteInches = regexp(tikzText, ...
    '(?m)^\s*(width|height)\s*=\s*[0-9.]+\s*in\s*,?', 'once');
assert(isempty(absoluteInches));

exportSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'exportPaperFigure.m'));
preflightSource = fileread(fullfile(projectRoot, 'toolbox', 'sweep', ...
    'preflightPaperFigureToolchain.m'));
configSource = fileread(fullfile(projectRoot, 'config', ...
    'getFairDOMPComparisonConfig.m'));
runnerSource = fileread(fullfile(projectRoot, 'run_parameter_sweep.m'));
selectedSource = fileread(fullfile(projectRoot, 'run_selected_comparison.m'));
matlabSources = string(exportSource) + newline + ...
    string(preflightSource) + newline + string(configSource) + newline + ...
    string(runnerSource) + newline + string(selectedSource);
assert(~contains(matlabSources, 'latexmk', 'IgnoreCase', true));
assert(~contains(matlabSources, 'resolveLatexmkCommand'));
assert(~contains(exportSource, 'system('));
assert(~contains(preflightSource, 'system('));
assert(~contains(exportSource, 'standaloneTex'));
assert(~contains(exportSource, 'standalonePdf'));
assert(~contains(exportSource, 'writeStandaloneWrapper'));
assert(~contains(exportSource, 'removeTikzLegendLayout'));
assert(~contains(exportSource, 'readlines('));
assert(~contains(exportSource, 'writelines('));

figureHandle = openfig(char(files.fig), 'new', 'visible');
drawnow;
assert(isgraphics(figureHandle, 'figure'));
assert(strcmpi(figureHandle.Visible, 'on'));
close(figureHandle);

clear cleanup;
fprintf('PAPER FIGURE ROUNDTRIP TEST: PASS\n');

function removeFixture(directory)
close all force;
if isfolder(directory)
    rmdir(directory, 's');
end
end
