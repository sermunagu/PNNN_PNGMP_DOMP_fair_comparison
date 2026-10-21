function files = exportPaperFigure(figureHandle, baseFilename, options)
% exportPaperFigure - Export one figure handle to FIG, PNG, and TikZ.
% MATLAB produces reusable source and preview artifacts. Final dimensions,
% fonts, legends, IEEEtran integration, and PDF compilation belong to Overleaf.

if nargin < 3 || isempty(options)
    options = struct();
end
baseFilename = char(string(baseFilename));
[directory, name] = fileparts(baseFilename);
if isempty(directory)
    directory = pwd;
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
visibilityCleanup = onCleanup(@() set(figureHandle, 'Visible', originalVisible));
figureHandle.Visible = 'on';
savefig(figureHandle, char(files.fig));
clear visibilityCleanup;
exportgraphics(figureHandle, char(files.png), 'Resolution', 300, ...
    'BackgroundColor', 'white');

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
