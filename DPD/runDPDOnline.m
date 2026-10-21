function [outputFile, dpd] = runDPDOnline(modelFile, inputFile)

% Load saved models, infer from x and save three DPD signals. No training.

root = fileparts(fileparts(mfilename('fullpath')));

addpath(genpath(fullfile(root, 'toolbox')));

M = load(modelFile, 'bundle');
S = load(inputFile, 'x');

dpd = applyDPDModels(M.bundle, S.x(:));

outDir = fileparts(modelFile);

if isempty(outDir), outDir = pwd; end

outputFile = [tempname(outDir) '.mat'];

save(outputFile, 'dpd', '-v7.3');
fprintf('Signals saved: %s\n', outputFile);

end
