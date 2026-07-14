function hComplex = coupledRealToComplexCoefficients(hReal, nComplex)
% coupledRealToComplexCoefficients - Recover complex GMP coefficients.
% The input layout is all real coefficient parts followed by all imaginary
% parts, matching the coupled real system used throughout PN-GMP phase 1.

if nargin < 2 || isempty(nComplex)
    nComplex = numel(hReal) / 2;
end
if ~isscalar(nComplex) || ~isfinite(nComplex) || ...
        nComplex < 1 || nComplex ~= floor(nComplex)
    error('coupledRealToComplexCoefficients:InvalidCount', ...
        'nComplex must be a positive integer scalar.');
end
if ~isnumeric(hReal) || ~isreal(hReal) || ~isvector(hReal) || ...
        numel(hReal) ~= 2*nComplex || any(~isfinite(hReal))
    error('coupledRealToComplexCoefficients:InvalidCoefficients', ...
        'hReal must be a finite real vector with exactly 2*nComplex values.');
end

hReal = hReal(:);
hComplex = complex(hReal(1:nComplex), hReal(nComplex+1:end));
end
