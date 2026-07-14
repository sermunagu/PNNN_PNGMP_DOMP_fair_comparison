function options = validateDOMPOptions(options)
% validateDOMPOptions - Validate the fixed DOMP numerical options.
% DOMP is the only sparse selector supported by the comparison projects;
% this helper provides tolerances without exposing an algorithm switch.

if nargin < 1 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isscalar(options)
    error('validateDOMPOptions:InvalidOptions', ...
        'options must be a scalar struct.');
end

defaults = struct( ...
    'columnTolerance', 1e-12, ...
    'correlationTolerance', 1e-14, ...
    'residualTolerance', 1e-12, ...
    'lsTolerance', []);
allowed = fieldnames(defaults);
unknown = setdiff(fieldnames(options), allowed);
if ~isempty(unknown)
    error('validateDOMPOptions:UnknownOption', ...
        'Unknown DOMP option: %s.', unknown{1});
end

for index = 1:numel(allowed)
    name = allowed{index};
    if ~isfield(options, name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end

validateNonnegativeScalar(options.columnTolerance, 'columnTolerance');
validateNonnegativeScalar(options.correlationTolerance, ...
    'correlationTolerance');
validateNonnegativeScalar(options.residualTolerance, 'residualTolerance');
if ~isempty(options.lsTolerance)
    validateNonnegativeScalar(options.lsTolerance, 'lsTolerance');
end
end

function validateNonnegativeScalar(value, name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value < 0
    error('validateDOMPOptions:InvalidTolerance', ...
        '%s must be a finite non-negative real scalar.', name);
end
end
