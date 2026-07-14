function U = buildGMPRegressorRows(x, rows, rManagerGMP, support)
% buildGMPRegressorRows - Evaluate an existing GMP population on selected rows.
% The helper preserves the periodic indexing and regressor ordering used by
% the current block GMP baseline without mutating the regressor manager.

x = x(:);
if isempty(x) || ~isnumeric(x) || ~isfloat(x) || any(~isfinite(x))
    error('buildGMPRegressorRows:InvalidInput', ...
        'x must be a non-empty finite floating-point vector.');
end

if nargin < 2 || isempty(rows)
    rows = (1:numel(x)).';
end
rows = validateIndexVector(rows, numel(x), 'rows');

if nargin < 3 || ~isobject(rManagerGMP) || ...
        ~isprop(rManagerGMP, 'regPopulation')
    error('buildGMPRegressorRows:InvalidManager', ...
        'rManagerGMP must expose a regPopulation property.');
end
regPopulation = rManagerGMP.regPopulation;
nPopulation = numel(regPopulation);
if nPopulation == 0
    error('buildGMPRegressorRows:EmptyPopulation', ...
        'rManagerGMP.regPopulation must not be empty.');
end

if nargin < 4 || isempty(support)
    support = 1:nPopulation;
end
support = validateIndexVector(support, nPopulation, 'support').';
if numel(unique(support, 'stable')) ~= numel(support)
    error('buildGMPRegressorRows:DuplicateSupport', ...
        'support must not contain duplicate indices.');
end

regSpecs = prepareRegressorSpecs(regPopulation);
shifts = collectShifts(regSpecs, support);
nRows = numel(rows);
nRegs = numel(support);
nSignal = numel(x);

tap = cell(numel(shifts), 1);
absTap = cell(numel(shifts), 1);
for k = 1:numel(shifts)
    idx = wrapIndex(rows - shifts(k), nSignal);
    tap{k} = x(idx);
    absTap{k} = abs(tap{k});
end

U = complex(zeros(nRows, nRegs));
for k = 1:nRegs
    spec = regSpecs(support(k));
    value = complex(ones(nRows, 1));

    for it = 1:numel(spec.Xq)
        tapValue = tap{find(shifts == spec.Xq(it), 1)};
        value = value .* (tapValue .^ spec.Xpow(it));
    end
    for it = 1:numel(spec.Xconjq)
        tapValue = tap{find(shifts == spec.Xconjq(it), 1)};
        value = value .* (conj(tapValue) .^ spec.Xconjpow(it));
    end
    for it = 1:numel(spec.Xenvq)
        envelopeValue = absTap{find(shifts == spec.Xenvq(it), 1)};
        value = value .* (envelopeValue .^ spec.Xenvpow(it));
    end

    U(:, k) = value;
end

if any(~isfinite(U), 'all')
    error('buildGMPRegressorRows:NonFiniteOutput', ...
        'The evaluated GMP matrix contains NaN or Inf values.');
end
end

function regSpecs = prepareRegressorSpecs(regPopulation)
nRegs = numel(regPopulation);
emptySpec = struct( ...
    'Xq', [], 'Xpow', [], ...
    'Xconjq', [], 'Xconjpow', [], ...
    'Xenvq', [], 'Xenvpow', []);
regSpecs = repmat(emptySpec, 1, nRegs);

for k = 1:nRegs
    [regSpecs(k).Xq, regSpecs(k).Xpow] = groupedTerms(regPopulation(k).X);
    [regSpecs(k).Xconjq, regSpecs(k).Xconjpow] = ...
        groupedTerms(regPopulation(k).Xconj);
    [regSpecs(k).Xenvq, regSpecs(k).Xenvpow] = ...
        groupedTerms(regPopulation(k).Xenv);
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

function shifts = collectShifts(regSpecs, support)
shifts = [];
for k = support(:).'
    shifts = [shifts, regSpecs(k).Xq, regSpecs(k).Xconjq, ...
        regSpecs(k).Xenvq]; %#ok<AGROW>
end

if isempty(shifts)
    shifts = 0;
else
    shifts = unique(shifts, 'stable');
end
end

function idx = wrapIndex(idx, nSignal)
idx = mod(idx - 1, nSignal) + 1;
end

function idx = validateIndexVector(idx, upperBound, name)
if ~isnumeric(idx) || ~isreal(idx) || ~isvector(idx) || isempty(idx) || ...
        any(~isfinite(idx)) || any(idx ~= floor(idx)) || ...
        any(idx < 1) || any(idx > upperBound)
    error('buildGMPRegressorRows:InvalidIndices', ...
        '%s must contain finite integer indices in [1, %d].', ...
        name, upperBound);
end
idx = double(idx(:));
end
