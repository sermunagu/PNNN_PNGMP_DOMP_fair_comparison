function selection = selectSharedIQFeatures( ...
    features, y_complex, maximum_features, domp_options)
% selectSharedIQFeatures - Select shared independent I/Q features by DOMP.
% A real feature matrix and complex phase-normalized target are passed to
% the central selector; I and Q coefficients are fitted after selection.

if nargin < 4
    domp_options = struct();
end
if ~isnumeric(features) || ~isreal(features) || ~isfloat(features) || ...
        isempty(features) || any(~isfinite(features), 'all')
    error('selectSharedIQFeatures:InvalidFeatures', ...
        'features must be a finite non-empty floating-point real matrix.');
end
y_complex = y_complex(:);
if numel(y_complex) ~= size(features, 1) || any(~isfinite(y_complex))
    error('selectSharedIQFeatures:InvalidTarget', ...
        'y_complex must contain one finite value per feature row.');
end
if ~isscalar(maximum_features) || maximum_features < 1 || ...
        maximum_features ~= floor(maximum_features) || ...
        maximum_features > size(features, 2)
    error('selectSharedIQFeatures:InvalidCount', ...
        'maximum_features must be a valid positive integer.');
end

[support, history] = selectDOMPSupport( ...
    features, y_complex, maximum_features, domp_options);
if numel(support) ~= maximum_features
    error('selectSharedIQFeatures:EarlyTermination', ...
        'DOMP selected %d of %d requested shared features.', ...
        numel(support), maximum_features);
end

selection = struct();
selection.supportFeatures = support(:);
selection.maximumFeatures = maximum_features;
selection.realParameters = 2*maximum_features;
selection.targetOutputs = 2;
selection.selectionDomainRows = size(features, 1);
selection.selectionMethod = 'DOMP';
selection.dompHistory = history;
end
