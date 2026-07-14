function yhat = GMP_blockPredict(x, rManagerGMP, support, h, blockSize, rows)
% GMP_blockPredict - Predict GMP output by row blocks.
%
% This helper evaluates the selected GMP regressors without materializing
% the full regression matrix, supporting memory-safe GMP baselines for PNNN.
%
% Inputs:
%   x - Full modeled-block input signal.
%   rManagerGMP, support, h - GMP basis manager, selected regressors, coefficients.
%   blockSize, rows - Block length and optional full-domain rows to evaluate.
%
% Outputs:
%   yhat - GMP prediction for the requested rows.

x = x(:);
N = numel(x);
if nargin < 6 || isempty(rows)
    rows = (1:N).';
else
    rows = rows(:);
end
if nargin < 5 || isempty(blockSize)
    blockSize = 8192;
end

support = support(:).';
h = h(:);
if numel(support) ~= numel(h)
    error('support y h deben tener la misma longitud.');
end
if isempty(support)
    error('support no puede estar vacio.');
end
if min(rows) < 1 || max(rows) > N
    error('rows fuera de rango para N=%d.', N);
end
if max(support) > numel(rManagerGMP.regPopulation)
    error('support contiene indices fuera de regPopulation.');
end

regSpecs = prepareRegressorSpecs(rManagerGMP.regPopulation);
blockSize = max(1, min(blockSize, numel(rows)));
yhat = complex(zeros(numel(rows), 1));

for first = 1:blockSize:numel(rows)
    last = min(first + blockSize - 1, numel(rows));
    blockRows = rows(first:last);
    Ublk = buildRegressorBlock(x, blockRows, regSpecs, support);
    yhat(first:last) = Ublk * h;
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
