function support = selectDOMPSupport(X, y, maxComponents, columnTolerance)
% selectDOMPSupport - Select a sparse support using fixed DOMP.
% Remaining candidates are orthogonalized with Gram-Schmidt; coefficients
% and the residual are recomputed from the original design matrix.

if isvector(y)
    y = y(:);
end

nCandidates = size(X, 2);
maxComponents = min(double(maxComponents), nCandidates);
Z = X;
residual = y;
selected = false(nCandidates, 1);
support = zeros(maxComponents, 1);

initialNorms = sqrt(sum(abs(X).^2, 1));
absoluteColumnTolerance = columnTolerance * max(1, max(initialNorms));
selectedCount = 0;

for iteration = 1:maxComponents
    columnNorms = sqrt(sum(abs(Z).^2, 1));
    eligible = ~selected.' & columnNorms > absoluteColumnTolerance;
    if ~any(eligible)
        break;
    end

    Z(:, eligible) = Z(:, eligible) ./ columnNorms(eligible);
    Z(:, ~eligible) = 0;

    correlations = Z(:, eligible)' * residual;
    scores = sqrt(sum(abs(correlations).^2, 2));
    eligibleIndices = find(eligible);
    [~, localIndex] = max(scores);
    bestIndex = eligibleIndices(localIndex);

    q = Z(:, bestIndex);
    Z = Z - q * (q' * Z);
    selected(bestIndex) = true;
    Z(:, bestIndex) = 0;

    selectedCount = selectedCount + 1;
    support(selectedCount) = bestIndex;
    activeSupport = support(1:selectedCount);
    coefficients = solveSelectedLeastSquares(X(:, activeSupport), y);
    residual = y - X(:, activeSupport) * coefficients;
end

support = support(1:selectedCount);
if isempty(support)
    error('selectDOMPSupport:EmptySupport', ...
        'DOMP did not find an eligible component.');
end
end




function coefficients = solveSelectedLeastSquares(XSelected, y)
% Preserve the baseline QR solve and minimum-norm rank-deficient fallback.
[Q, R] = qr(XSelected, 0);
tolerance = max(size(XSelected))*eps(norm(R, 2));
diagonal = abs(diag(R));
if isempty(diagonal) || any(diagonal <= tolerance)
    coefficients = lsqminnorm(XSelected, y);
else
    coefficients = R \ (Q' * y);
end
end
