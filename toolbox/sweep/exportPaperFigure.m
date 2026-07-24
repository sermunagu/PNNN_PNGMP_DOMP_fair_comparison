function files = exportPaperFigure(figureHandle, baseFilename, options)
% exportPaperFigure - Export one figure handle to FIG, PNG, and TikZ.
% MATLAB produces reusable source and preview artifacts. Final dimensions,
% fonts, legends, IEEEtran integration, and PDF compilation belong to Overleaf.

if nargin < 3 || isempty(options)
    options = struct();
end
if ~isgraphics(figureHandle, 'figure')
    error('exportPaperFigure:InvalidFigure', ...
        'figureHandle must refer to one live MATLAB figure.');
end
baseFilename = char(string(baseFilename));
[directory, name, extension] = fileparts(baseFilename);
if isempty(directory)
    directory = pwd;
end
if ~isempty(extension)
    error('exportPaperFigure:BaseFilenameHasExtension', ...
        'baseFilename must not include a file extension.');
end
if ~isfolder(directory)
    mkdir(directory);
end
files = struct( ...
    'fig', string(fullfile(directory, name + ".fig")), ...
    'png', string(fullfile(directory, name + ".png")), ...
    'tikz', string(fullfile(directory, name + ".tikz")));

set(figureHandle, 'Color', 'w', 'InvertHardcopy', 'off');
originalVisible = figureHandle.Visible;
visibilityCleanup = onCleanup(@() restoreVisibility( ...
    figureHandle, originalVisible));
figureHandle.Visible = 'on';
savefig(figureHandle, char(files.fig));
clear visibilityCleanup;
exportgraphics(figureHandle, char(files.png), 'Resolution', 300, ...
    'BackgroundColor', 'white');

if exist('matlab2tikz', 'file') ~= 2
    error('exportPaperFigure:MissingMatlab2tikz', ...
        ['matlab2tikz is not on the MATLAB path. Initialize the submodule ' ...
        'and add third_party/matlab2tikz/src.']);
end
matlab2tikzArguments = {char(files.tikz), ...
    'figurehandle', figureHandle, ...
    'width', '\figurewidth', 'height', '\figureheight', ...
    'standalone', false, 'showInfo', false};
if isfield(options, 'tikzExtraAxisOptions')
    matlab2tikzArguments = [matlab2tikzArguments, ...
        {'extraAxisOptions', options.tikzExtraAxisOptions}];
end
matlab2tikz(matlab2tikzArguments{:});
end

function restoreVisibility(figureHandle, originalVisible)
if isgraphics(figureHandle, 'figure')
    figureHandle.Visible = originalVisible;
end
end
