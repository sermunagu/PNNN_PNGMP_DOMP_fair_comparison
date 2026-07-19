function report = preflightPaperFigureToolchain(projectRoot, paperConfig)
% preflightPaperFigureToolchain - Fail before data loading or training.
% A standalone pgfplots document is compiled with the same LaTeX packages
% used by exportPaperFigure so missing dependencies are detected up front.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
if nargin < 2 || isempty(paperConfig)
    paperConfig = struct();
end
if ~isfield(paperConfig, 'matlab2tikzSource')
    paperConfig.matlab2tikzSource = fullfile(projectRoot, ...
        'third_party', 'matlab2tikz', 'src');
end
if ~isfield(paperConfig, 'latexmkCommand')
    paperConfig.latexmkCommand = 'latexmk';
end

sourceDirectory = char(string(paperConfig.matlab2tikzSource));
if ~isfile(fullfile(sourceDirectory, 'matlab2tikz.m'))
    error('preflightPaperFigureToolchain:MissingSubmodule', ...
        ['matlab2tikz was not found at %s. Run: git submodule update ' ...
        '--init --recursive'], sourceDirectory);
end
addpath(sourceDirectory);
if exist('matlab2tikz', 'file') ~= 2
    error('preflightPaperFigureToolchain:MissingMatlab2tikz', ...
        'matlab2tikz is not callable after adding %s.', sourceDirectory);
end

latexmkCommand = char(string(paperConfig.latexmkCommand));
try
    latexmk = resolveLatexmkCommand(latexmkCommand);
catch exception
    error('preflightPaperFigureToolchain:MissingLatexmk', ...
        'latexmk is unavailable: %s', exception.message);
end

compileDirectory = tempname;
mkdir(compileDirectory);
tikzFile = fullfile(compileDirectory, 'paper_figure_preflight.tikz');
texFile = fullfile(compileDirectory, 'paper_figure_preflight.tex');
figureHandle = figure('Visible', 'off', 'Color', 'w');
figureCleanup = onCleanup(@() close(figureHandle));
style = getIEEEPaperStyle();
plot(axes(figureHandle), [0 1], [0 1], ...
    'Color', style.gmpBlue, 'LineWidth', style.mainLineWidth);
matlab2tikz(tikzFile, 'figurehandle', figureHandle, ...
    'width', '\figurewidth', 'height', '\figureheight', ...
    'standalone', false, 'showInfo', false);
clear figureCleanup;
writePreflightDocument(texFile, 'paper_figure_preflight.tikz');
command = sprintf(['%s -pdf -interaction=nonstopmode ' ...
    '-halt-on-error -outdir="%s" "%s"'], latexmk.commandPrefix, ...
    compileDirectory, texFile);
[status, compileText] = system(command);
logFile = fullfile(compileDirectory, 'paper_figure_preflight.log');
pdfFile = fullfile(compileDirectory, 'paper_figure_preflight.pdf');
if status ~= 0 || ~isfile(pdfFile)
    error('preflightPaperFigureToolchain:LatexCompileFailed', ...
        ['The minimal standalone pgfplots compile failed before training. ' ...
        'Inspect %s. latexmk output: %s'], logFile, strtrim(compileText));
end
rmdir(compileDirectory, 's');

versionLines = splitlines(strtrim(latexmk.versionText));
versionRow = find(contains(versionLines, 'Latexmk,'), 1, 'first');
if isempty(versionRow)
    versionRow = find(strlength(versionLines) > 0, 1, 'first');
end
firstLine = versionLines(versionRow);
report = struct('matlab2tikzSource', string(sourceDirectory), ...
    'latexmkCommand', string(latexmkCommand), ...
    'latexmkCommandPrefix', latexmk.commandPrefix, ...
    'usedMatlabPerlFallback', latexmk.usedMatlabPerlFallback, ...
    'latexmkVersion', firstLine, 'minimalCompilePassed', true);
end

function writePreflightDocument(filename, tikzFilename)
lines = [ ...
    "\documentclass[tikz,border=2pt]{standalone}"; ...
    "\usepackage{pgfplots}"; ...
    "\usepackage{tikz}"; ...
    "\usepackage{amsmath}"; ...
    "\pgfplotsset{compat=1.18}"; ...
    "\newlength\figurewidth"; ...
    "\newlength\figureheight"; ...
    "\setlength\figurewidth{3.5in}"; ...
    "\setlength\figureheight{2.4in}"; ...
    "\begin{document}"; ...
    "\input{" + string(tikzFilename) + "}"; ...
    "\end{document}"];
writelines(lines, filename);
end
