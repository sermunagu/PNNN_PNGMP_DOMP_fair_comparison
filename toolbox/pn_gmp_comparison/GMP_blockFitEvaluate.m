function [res, fit] = GMP_blockFitEvaluate( ...
    x, y, identificationIndices, validationIndices, ...
    rManagerGMP, cfg, label)
% GMP_blockFitEvaluate - Fit and evaluate a GMP basis by row blocks.
%
% The helper accumulates normal equations only over identification rows and
% evaluates the fitted models over a separately supplied validation domain.
%
% Inputs:
%   x, y - Full modeled-block input and output signals.
%   identificationIndices - Rows used for DOMP and coefficient fitting.
%   validationIndices - Rows used for full-signal validation.
%   rManagerGMP, cfg, label - GMP basis manager, config, and log label.
%
% Outputs:
%   res - Struct with pinv/ridge NMSE metrics and diagnostics.
%   fit - Struct with selected support and coefficient fits.

if nargin < 7 || isempty(label)
    label = 'GMP block';
end
if nargin < 6 || isempty(cfg)
    cfg = struct();
end

x = x(:);
y = y(:);
N = numel(x);
if isempty(x) || numel(y) ~= N
    error('x e y deben ser vectores no vacios de la misma longitud.');
end
if any(~isfinite(x)) || any(~isfinite(y))
    error('x e y deben contener valores finitos.');
end

identificationIndices = identificationIndices(:);
validationIndices = validationIndices(:);
if isempty(identificationIndices) || isempty(validationIndices)
    error('Identification and validation indices must not be empty.');
end
if min([identificationIndices; validationIndices]) < 1 || ...
        max([identificationIndices; validationIndices]) > N
    error('GMP indices are outside the signal domain: N=%d min=%d max=%d.', ...
        N, min([identificationIndices; validationIndices]), ...
        max([identificationIndices; validationIndices]));
end

blockSize = getCfgField(cfg, 'blockSize', 8192);
selectionMethod = string(getCfgField(cfg, 'selectionMethod', 'DOMP'));
if selectionMethod ~= "DOMP"
    error('GMP_blockFitEvaluate:InvalidSelectionMethod', ...
        'DOMP is the only supported sparse selection method.');
end
dompOptions = getCfgField(cfg, 'dompOptions', struct());
maxPopulation = getCfgField(cfg, 'maxPopulation', 100);
lambda1 = getCfgField(cfg, 'lambda1', 1e-3);
lambda2 = getCfgField(cfg, 'lambda2', 1e-4);

regSpecs = prepareRegressorSpecs(rManagerGMP.regPopulation);
nRegsTotal = numel(regSpecs);
if nRegsTotal == 0
    error('La poblacion GMP esta vacia.');
end

blockSize = max(1, min(blockSize, ...
    max(numel(identificationIndices), numel(validationIndices))));
maxPopulation = max(1, min(maxPopulation, nRegsTotal));

bytesComplex = 16;
fullUGB = N * nRegsTotal * bytesComplex / 2^30;
oldMatrixGB = 3 * fullUGB; % Historical full-matrix materialization estimate.
blockUGB = blockSize * nRegsTotal * bytesComplex / 2^30;
gramGB = nRegsTotal * nRegsTotal * bytesComplex / 2^30;
selectionUGB = numel(identificationIndices) * nRegsTotal * ...
    bytesComplex / 2^30;
newPeakGB = 2*selectionUGB + blockUGB + gramGB;

fprintf(['[%s] N=%d | identification=%d | full validation=%d | ' ...
    'regresores=%d | blockSize=%d\n'], ...
    label, N, numel(identificationIndices), numel(validationIndices), ...
    nRegsTotal, blockSize);
fprintf('[%s] U completa estimada: %.2f GB (NO se materializa)\n', label, fullUGB);
fprintf('[%s] Pico viejo estimado solo matrices U/Un: %.2f GB\n', label, oldMatrixGB);
fprintf(['[%s] U_identification DOMP: %.3f GB | U_block: %.3f GB | ' ...
    'pico analitico aprox: %.3f GB\n'], ...
    label, selectionUGB, blockUGB, newPeakGB);

tStart = tic;
fprintf('[%s] Construyendo la matriz de identificacion para DOMP...\n', label);
UIdentification = buildDesignMatrix( ...
    x, identificationIndices, regSpecs, 1:nRegsTotal, blockSize);
yIdentification = y(identificationIndices);
[support, dompHistory] = selectDOMPSupport( ...
    UIdentification, yIdentification, maxPopulation, dompOptions);
normU = sqrt(sum(abs(UIdentification).^2, 1)).';
normU(normU == 0) = 1;
USelectedNormalized = UIdentification(:, support) ./ normU(support).';
Gs = USelectedNormalized' * USelectedNormalized;
bs = USelectedNormalized' * yIdentification;
y2Identification = sum(abs(yIdentification).^2);
nmsePath = 20*log10(max(dompHistory.residualNorm, realmin) ./ ...
    max(sqrt(y2Identification), realmin));
clear UIdentification USelectedNormalized yIdentification;

nActive = numel(support);
fprintf('[%s] Regresores activos: %d/%d (selectionMethod=%s)\n', ...
    label, nActive, nRegsTotal, selectionMethod);
I = eye(nActive);

aLS = Gs \ bs;
aRidge1 = (Gs + lambda1 * I) \ bs;
aRidge2 = (Gs + lambda2 * I) \ bs;

hLS = aLS ./ normU(support);
hRidge1 = aRidge1 ./ normU(support);
hRidge2 = aRidge2 ./ normU(support);

fprintf('[%s] Evaluando NMSE por bloques...\n', label);
H = [hLS(:), hRidge1(:), hRidge2(:)];
[err2IdentificationAll, y2IdentificationEval] = ...
    predictionErrorEnergyMulti(x, y, identificationIndices, ...
    regSpecs, support, H, blockSize);
[err2ValidationAll, y2Validation] = ...
    predictionErrorEnergyMulti(x, y, validationIndices, ...
    regSpecs, support, H, blockSize);
err2IdentificationLS = err2IdentificationAll(1);
err2IdentificationR1 = err2IdentificationAll(2);
err2IdentificationR2 = err2IdentificationAll(3);
err2ValidationLS = err2ValidationAll(1);
err2ValidationR1 = err2ValidationAll(2);
err2ValidationR2 = err2ValidationAll(3);

elapsed = toc(tStart);

res = struct();
res.NMSE_identification_pinv = nmseFromEnergy( ...
    err2IdentificationLS, y2IdentificationEval);
res.NMSE_validation_full_signal_pinv = nmseFromEnergy( ...
    err2ValidationLS, y2Validation);
res.NMSE_identification_ridge_1e3 = nmseFromEnergy( ...
    err2IdentificationR1, y2IdentificationEval);
res.NMSE_validation_full_signal_ridge_1e3 = nmseFromEnergy( ...
    err2ValidationR1, y2Validation);
res.NMSE_identification_ridge_1e4 = nmseFromEnergy( ...
    err2IdentificationR2, y2IdentificationEval);
res.NMSE_validation_full_signal_ridge_1e4 = nmseFromEnergy( ...
    err2ValidationR2, y2Validation);
res.N_signal = N;
res.N = N;
res.N_identification = numel(identificationIndices);
res.N_validation = numel(validationIndices);
res.Un_rows = N;
res.L = N;
res.Un_cols = nRegsTotal;
res.nRegressorsTotal = nRegsTotal;
res.nCoeff_GMP = nActive;
res.blockSize = blockSize;
res.selectionMethod = char(selectionMethod);
res.support = support(:);
res.dompHistory = dompHistory;
res.elapsedSeconds = elapsed;
res.estimatedFullUGB = fullUGB;
res.estimatedOldMatrixGB = oldMatrixGB;
res.estimatedBlockUGB = blockUGB;
res.estimatedGramGB = gramGB;
res.estimatedNewPeakGB = newPeakGB;
res.maxCoeff_pinv = max(abs(hLS));
res.maxCoeff_ridge_1e3 = max(abs(hRidge1));
res.maxCoeff_ridge_1e4 = max(abs(hRidge2));

fit = struct();
fit.support = support(:);
fit.nmsePath = nmsePath(:);
fit.h_pinv = hLS(:);
fit.h_ridge_1e3 = hRidge1(:);
fit.h_ridge_1e4 = hRidge2(:);
fit.normU = normU(:);
fit.selectionMethod = char(selectionMethod);
fit.dompHistory = dompHistory;

fprintf(['[%s] Terminado en %.2f s | NMSE full validation LS=%.2f dB | ' ...
    'ridge1e-3=%.2f dB | ridge1e-4=%.2f dB\n'], ...
    label, elapsed, res.NMSE_validation_full_signal_pinv, ...
    res.NMSE_validation_full_signal_ridge_1e3, ...
    res.NMSE_validation_full_signal_ridge_1e4);
end

function U = buildDesignMatrix(x, rows, regSpecs, regIdx, blockSize)
U = complex(zeros(numel(rows), numel(regIdx)));
for first = 1:blockSize:numel(rows)
    last = min(first + blockSize - 1, numel(rows));
    blockRows = rows(first:last);
    U(first:last, :) = buildRegressorBlock( ...
        x, blockRows, regSpecs, regIdx);
end
end

function [err2, y2] = predictionErrorEnergyMulti(x, y, rows, regSpecs, regIdx, H, blockSize)
err2 = zeros(1, size(H, 2));
y2 = 0;
for first = 1:blockSize:numel(rows)
    last = min(first + blockSize - 1, numel(rows));
    blockRows = rows(first:last);
    Ublk = buildRegressorBlock(x, blockRows, regSpecs, regIdx);
    yblk = y(blockRows);
    E = Ublk * H - yblk;
    err2 = err2 + sum(abs(E).^2, 1);
    y2 = y2 + sum(abs(yblk).^2);
end
end

function Ublk = buildRegressorBlock(x, rows, regSpecs, regIdx)
rows = rows(:);
nRows = numel(rows);
nRegs = numel(regIdx);
N = numel(x);

shifts = collectShifts(regSpecs, regIdx);
tap = cell(numel(shifts), 1);
abstap = cell(numel(shifts), 1);
for k = 1:numel(shifts)
    idx = wrapIndex(rows - shifts(k), N);
    tap{k} = x(idx);
    abstap{k} = abs(tap{k});
end

Ublk = complex(zeros(nRows, nRegs));
for k = 1:nRegs
    spec = regSpecs(regIdx(k));
    v = complex(ones(nRows, 1));
    for it = 1:numel(spec.Xq)
        tv = tap{find(shifts == spec.Xq(it), 1)};
        v = v .* (tv .^ spec.Xpow(it));
    end
    for it = 1:numel(spec.Xconjq)
        tv = tap{find(shifts == spec.Xconjq(it), 1)};
        v = v .* (conj(tv) .^ spec.Xconjpow(it));
    end
    for it = 1:numel(spec.Xenvq)
        av = abstap{find(shifts == spec.Xenvq(it), 1)};
        v = v .* (av .^ spec.Xenvpow(it));
    end
    Ublk(:, k) = v;
end
end

function regSpecs = prepareRegressorSpecs(regPopulation)
nRegs = numel(regPopulation);
emptySpec = struct('Xq', [], 'Xpow', [], 'Xconjq', [], 'Xconjpow', [], 'Xenvq', [], 'Xenvpow', []);
regSpecs = repmat(emptySpec, 1, nRegs);
for k = 1:nRegs
    [regSpecs(k).Xq, regSpecs(k).Xpow] = groupedTerms(regPopulation(k).X);
    [regSpecs(k).Xconjq, regSpecs(k).Xconjpow] = groupedTerms(regPopulation(k).Xconj);
    [regSpecs(k).Xenvq, regSpecs(k).Xenvpow] = groupedTerms(regPopulation(k).Xenv);
end
end

function [q, p] = groupedTerms(terms)
terms = terms(:).';
if isempty(terms)
    q = [];
    p = [];
    return;
end
q = unique(terms, 'stable');
p = zeros(size(q));
for k = 1:numel(q)
    p(k) = sum(terms == q(k));
end
end

function shifts = collectShifts(regSpecs, regIdx)
shifts = [];
for k = regIdx(:).'
    shifts = [shifts, regSpecs(k).Xq, regSpecs(k).Xconjq, regSpecs(k).Xenvq]; %#ok<AGROW>
end
if isempty(shifts)
    shifts = 0;
else
    shifts = unique(shifts, 'stable');
end
end

function idx = wrapIndex(idx, N)
idx = mod(idx - 1, N) + 1;
end

function nmse = nmseFromEnergy(err2, y2)
nmse = 10 * log10(max(real(err2), realmin) / max(real(y2), realmin));
end

function value = getCfgField(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
