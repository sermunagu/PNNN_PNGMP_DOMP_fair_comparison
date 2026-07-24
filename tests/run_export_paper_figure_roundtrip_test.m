% Verify that exportPaperFigure writes a FIG readable by openfig.

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

assert(isfile(files.tikz));
assert(isfile(files.standaloneTex));
tikzText = fileread(files.tikz);
wrapperText = fileread(files.standaloneTex);
assert(contains(tikzText, 'width=0.68\figurewidth'));
assert(contains(tikzText, 'height=0.7\figureheight'));
assert(contains(wrapperText, ...
    '\PassOptionsToClass{lettersize,journal}{IEEEtran}'));

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
