function report = preflightPaperFigureToolchain(projectRoot, paperConfig)
% preflightPaperFigureToolchain - Verify matlab2tikz export only.
% LaTeX classes, wrappers, and PDF compilation are intentionally delegated to
% Overleaf and are not prerequisites for running the MATLAB pipeline.

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

exportDirectory = tempname;
mkdir(exportDirectory);
directoryCleanup = onCleanup(@() removeDirectory(exportDirectory));
tikzFile = fullfile(exportDirectory, 'paper_figure_preflight.tikz');
figureHandle = figure('Visible', 'off', 'Color', 'w');
figureCleanup = onCleanup(@() close(figureHandle));
style = getIEEEPaperStyle();
plot(axes(figureHandle), [0 1], [0 1], ...
    'Color', style.gmpBlue, 'LineWidth', style.mainLineWidth);
matlab2tikz(tikzFile, 'figurehandle', figureHandle, ...
    'width', '\figurewidth', 'height', '\figureheight', ...
    'standalone', false, 'showInfo', false);
if ~isfile(tikzFile) || dir(tikzFile).bytes == 0
    error('preflightPaperFigureToolchain:TikzExportFailed', ...
        'matlab2tikz did not produce a nonempty TikZ file.');
end

report = struct('matlab2tikzSource', string(sourceDirectory), ...
    'tikzExportPassed', true);
clear figureCleanup directoryCleanup;
end

function removeDirectory(directory)
if isfolder(directory)
    rmdir(directory, 's');
end
end
