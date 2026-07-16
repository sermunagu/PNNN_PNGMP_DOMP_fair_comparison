function savedConfig = loadOptionalExperimentSignature(filename)
% loadOptionalExperimentSignature - Read a signature only when it exists.
% Legacy comparison MAT files without the variable are ignored silently.
% The returned structure is suitable for strict reuse validation.

validPath = (ischar(filename) && isrow(filename)) || ...
    (isstring(filename) && isscalar(filename));
if ~validPath || ismissing(string(filename)) || ...
        strlength(string(filename)) == 0
    error('loadOptionalExperimentSignature:InvalidPath', ...
        'The MAT filename must be a nonempty text scalar.');
end
filename = char(string(filename));
savedConfig = struct();
if ~isfile(filename)
    return;
end
variables = whos('-file', filename);
if any(strcmp({variables.name}, 'experiment_signature'))
    savedConfig = load(filename, 'experiment_signature');
end
end
