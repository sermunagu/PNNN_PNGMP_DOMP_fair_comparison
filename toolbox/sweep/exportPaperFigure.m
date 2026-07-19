function files = exportPaperFigure(figureHandle, baseFilename, options)
% exportPaperFigure - Export one figure handle to FIG, PNG, TikZ, and PDF.
% The PDF is compiled from the same matlab2tikz output through a standalone
% wrapper. LaTeX intermediates are deleted only after successful compilation.

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
if ~isfield(options, 'latexmkCommand')
    options.latexmkCommand = 'latexmk';
end
if ~isfield(options, 'figureWidth')
    options.figureWidth = '7.0in';
end
if ~isfield(options, 'figureHeight')
    options.figureHeight = '4.5in';
end

files = struct( ...
    'fig', string(fullfile(directory, name + ".fig")), ...
    'png', string(fullfile(directory, name + ".png")), ...
    'tikz', string(fullfile(directory, name + ".tikz")), ...
    'pdf', string(fullfile(directory, name + ".pdf")));

set(figureHandle, 'Color', 'w', 'InvertHardcopy', 'off');
savefig(figureHandle, char(files.fig));
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

wrapperBase = name + "_standalone_wrapper";
wrapperTex = fullfile(directory, wrapperBase + ".tex");
writeStandaloneWrapper(wrapperTex, name + ".tikz", options);
latexmkCommand = char(string(options.latexmkCommand));
latexmk = resolveLatexmkCommand(latexmkCommand);
command = sprintf(['%s -pdf -interaction=nonstopmode ' ...
    '-halt-on-error -outdir="%s" "%s"'], latexmk.commandPrefix, ...
    directory, wrapperTex);
[status, commandOutput] = system(command);
wrapperPdf = fullfile(directory, wrapperBase + ".pdf");
wrapperLog = fullfile(directory, wrapperBase + ".log");
if status ~= 0 || ~isfile(wrapperPdf)
    error('exportPaperFigure:LatexCompileFailed', ...
        ['latexmk could not compile the standalone paper figure. ' ...
        'Intermediates were preserved; inspect %s. Output: %s'], ...
        wrapperLog, strtrim(commandOutput));
end

[moved, message] = movefile(wrapperPdf, char(files.pdf), 'f');
if ~moved
    error('exportPaperFigure:PDFMoveFailed', ...
        ['The compiled wrapper was preserved because the final PDF could ' ...
        'not be installed: %s'], message);
end
cleanupSuccessfulCompile(directory, wrapperBase);
end

function writeStandaloneWrapper(filename, tikzFilename, options)
lines = [ ...
    "\documentclass[tikz,border=2pt]{standalone}"; ...
    "\usepackage{pgfplots}"; ...
    "\usepackage{tikz}"; ...
    "\usepackage{amsmath}"; ...
    "\usetikzlibrary{plotmarks}"; ...
    "\usepgfplotslibrary{groupplots}"; ...
    "\pgfplotsset{compat=1.18}"; ...
    "\newlength\figurewidth"; ...
    "\newlength\figureheight"; ...
    "\setlength\figurewidth{" + string(options.figureWidth) + "}"; ...
    "\setlength\figureheight{" + string(options.figureHeight) + "}"; ...
    "\begin{document}"; ...
    "\input{" + replace(string(tikzFilename), "\", "/") + "}"; ...
    "\end{document}"];
writelines(lines, filename);
end

function cleanupSuccessfulCompile(directory, wrapperBase)
extensions = [".aux", ".fdb_latexmk", ".fls", ".log", ".tex"];
for extension = extensions
    filename = fullfile(directory, wrapperBase + extension);
    if isfile(filename)
        delete(filename);
    end
end
end
