function fit = fitPhaseNormGMPReal( ...
    U_phase_normalized, y_phase_normalized, fitCfg)
% fitPhaseNormGMPReal - Fit the coupled real form of a phase-normalized GMP.
% Phase 1 supports the full complex population or a frozen complex support,
% using the same paired column scaling as the current complex GMP baseline.

if nargin < 3 || isempty(fitCfg)
    fitCfg = struct();
end
if ~isstruct(fitCfg)
    error('fitPhaseNormGMPReal:InvalidConfig', ...
        'fitCfg must be a struct.');
end
if ~isnumeric(U_phase_normalized) || ...
        ~ismatrix(U_phase_normalized) || ...
        isempty(U_phase_normalized) || ...
        any(~isfinite(U_phase_normalized), 'all')
    error('fitPhaseNormGMPReal:InvalidMatrix', ...
        'U_phase_normalized must be a non-empty finite numeric matrix.');
end

nRows = size(U_phase_normalized, 1);
nComplexTotal = size(U_phase_normalized, 2);
if ~isnumeric(y_phase_normalized) || ...
        ~isvector(y_phase_normalized) || ...
        numel(y_phase_normalized) ~= nRows || ...
        any(~isfinite(y_phase_normalized))
    error('fitPhaseNormGMPReal:InvalidTarget', ...
        ['y_phase_normalized must be a finite vector with one value ' ...
         'per matrix row.']);
end
y_phase_normalized = y_phase_normalized(:);

supportMode = lower(string(getCfgField(fitCfg, 'supportMode', 'all')));
switch supportMode
    case "all"
        supportComplex = (1:nComplexTotal).';
    case "reuse_complex_support"
        supportComplex = getCfgField(fitCfg, 'supportComplex', []);
        supportComplex = validateSupport(supportComplex, nComplexTotal);
    otherwise
        error('fitPhaseNormGMPReal:InvalidSupportMode', ...
            'supportMode must be ''all'' or ''reuse_complex_support''.');
end

lambda = getCfgField(fitCfg, 'lambda', 0);
if ~isscalar(lambda) || ~isnumeric(lambda) || ~isreal(lambda) || ...
        ~isfinite(lambda) || lambda < 0
    error('fitPhaseNormGMPReal:InvalidLambda', ...
        'lambda must be a finite non-negative real scalar.');
end

if isfield(fitCfg, 'normUComplex') && ~isempty(fitCfg.normUComplex)
    normUComplex = fitCfg.normUComplex(:);
    if ~isnumeric(normUComplex) || ~isreal(normUComplex) || ...
            numel(normUComplex) ~= nComplexTotal || ...
            any(~isfinite(normUComplex)) || any(normUComplex <= 0)
        error('fitPhaseNormGMPReal:InvalidNorms', ...
            'normUComplex must contain one positive finite norm per complex column.');
    end
else
    normUComplex = sqrt(sum(abs(U_phase_normalized).^2, 1)).';
    normUComplex(normUComplex == 0) = 1;
end

U_selected = U_phase_normalized(:, supportComplex);
normSelected = normUComplex(supportComplex);
[U_real, y_real] = complexToCoupledReal( ...
    U_selected, y_phase_normalized);

normUReal = [normSelected; normSelected];
U_real_normalized = U_real ./ normUReal.';
GramRealNormalized = U_real_normalized.' * U_real_normalized;
rhsRealNormalized = U_real_normalized.' * y_real;

nActiveComplex = numel(supportComplex);
regularizedGram = GramRealNormalized + lambda * eye(2*nActiveComplex);
hNormalizedReal = regularizedGram \ rhsRealNormalized;
hReal = hNormalizedReal ./ normUReal;
hComplexRecovered = coupledRealToComplexCoefficients( ...
    hReal, nActiveComplex);

U_complex_normalized = U_selected ./ normSelected.';
GramComplexNormalized = U_complex_normalized' * U_complex_normalized;
if isfield(fitCfg, 'conditionNumber') && ~isempty(fitCfg.conditionNumber)
    conditionNumber = fitCfg.conditionNumber;
    if ~isscalar(conditionNumber) || ~isnumeric(conditionNumber) || ...
            ~isreal(conditionNumber) || isnan(conditionNumber) || ...
            conditionNumber < 0
        error('fitPhaseNormGMPReal:InvalidConditionNumber', ...
            'conditionNumber must be a non-negative real scalar.');
    end
else
    conditionNumber = cond(GramComplexNormalized);
end

supportReal = [supportComplex; nComplexTotal + supportComplex];
supportRealLocal = [(1:nActiveComplex).'; ...
    nActiveComplex + (1:nActiveComplex).'];
if numel(supportReal) ~= 2*nActiveComplex
    error('fitPhaseNormGMPReal:InvalidSupportPairing', ...
        'Every complex regressor must map to exactly two real columns.');
end

fit = struct();
fit.supportMode = char(supportMode);
fit.supportComplex = supportComplex;
fit.supportReal = supportReal;
fit.supportRealLocal = supportRealLocal;
fit.hReal = hReal;
fit.hComplexRecovered = hComplexRecovered;
fit.hComplexFull = zeros(nComplexTotal, 1, ...
    'like', hComplexRecovered);
fit.hComplexFull(supportComplex) = hComplexRecovered;
fit.hNormalizedReal = hNormalizedReal;
fit.normUComplex = normUComplex;
fit.normUReal = normUReal;
fit.lambda = double(lambda);
fit.solver = solverName(lambda);
fit.conditionNumber = conditionNumber;
fit.conditionNumberGram = conditionNumber;
fit.nComplexTotal = nComplexTotal;
fit.nActiveComplex = nActiveComplex;
fit.nActiveReal = 2*nActiveComplex;
if lambda == 0 && conditionNumber > 1e10
    fit.status = 'INCONCLUSIVE_CONDITIONING';
else
    fit.status = 'OK';
end
end

function support = validateSupport(support, upperBound)
if ~isnumeric(support) || ~isreal(support) || ~isvector(support) || ...
        isempty(support) || any(~isfinite(support)) || ...
        any(support ~= floor(support)) || ...
        any(support < 1) || any(support > upperBound)
    error('fitPhaseNormGMPReal:InvalidSupport', ...
        'supportComplex must contain integer indices in [1, %d].', upperBound);
end
support = double(support(:));
if numel(unique(support, 'stable')) ~= numel(support)
    error('fitPhaseNormGMPReal:DuplicateSupport', ...
        'supportComplex must not contain duplicate indices.');
end
end

function name = solverName(lambda)
if lambda == 0
    name = 'LS';
elseif lambda == 1e-3
    name = 'ridge_1e-3';
elseif lambda == 1e-4
    name = 'ridge_1e-4';
else
    name = sprintf('ridge_%g', lambda);
end
end

function value = getCfgField(cfg, fieldName, defaultValue)
if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end
end
