function [payload, reusable, filename] = updateSweepCheckpoint( ...
    resultDirectory, artifactKind, target, sweepIdentity, ...
    experimentSignature, newPayload)
% updateSweepCheckpoint - Load or atomically replace one sweep artifact.
% Artifact names are deterministic. Each call reads only the requested unit,
% so resuming one PNNN target never scans or rewrites previous targets.

if nargin < 6
    newPayload = [];
end
resultDirectory = normalizePath(resultDirectory, 'resultDirectory');
resultDirectoryPath = char(resultDirectory);
if ~isstruct(sweepIdentity) || ~isfield(sweepIdentity, 'digest') || ...
        ~validExperimentSignature(experimentSignature)
    error('updateSweepCheckpoint:InvalidIdentity', ...
        'Signed experiment and sweep identities are required.');
end
if ~isfolder(resultDirectoryPath)
    mkdir(resultDirectoryPath);
end
filename = string(fullfile( ...
    resultDirectoryPath, artifactName(artifactKind, target)));
filenamePath = char(filename);

if isempty(newPayload)
    [payload, reusable] = loadArtifact( ...
        filenamePath, sweepIdentity, experimentSignature);
    return;
end

checkpointArtifact = struct('schemaVersion', 1, ...
    'sweepIdentity', sweepIdentity, ...
    'experimentSignature', experimentSignature, ...
    'payload', newPayload);
[directory, ~, ~] = fileparts(filenamePath);
temporary = string(tempname(directory)) + ".mat";
temporaryPath = char(temporary);
cleanup = onCleanup(@() deleteIfPresent(temporaryPath));
save(temporaryPath, 'checkpointArtifact', '-v7.3');
[verified, reusable] = loadArtifact( ...
    temporaryPath, sweepIdentity, experimentSignature);
if ~reusable || isempty(verified)
    error('updateSweepCheckpoint:SerializationFailed', ...
        'The temporary sweep artifact failed verification.');
end
[moved, message] = movefile(temporaryPath, filenamePath, 'f');
if ~moved
    error('updateSweepCheckpoint:AtomicMoveFailed', ...
        'Could not install the verified sweep artifact: %s', message);
end
payload = newPayload;
reusable = true;
clear cleanup;
end

function [payload, reusable] = loadArtifact( ...
    filename, expectedIdentity, expectedSignature)
filename = normalizePath(filename, 'filename');
filenamePath = char(filename);
payload = [];
reusable = false;
if ~isfile(filenamePath)
    return;
end
try
    saved = load(filenamePath, 'checkpointArtifact');
    if ~isfield(saved, 'checkpointArtifact')
        return;
    end
    artifact = saved.checkpointArtifact;
    required = {'schemaVersion','sweepIdentity', ...
        'experimentSignature','payload'};
    if ~isstruct(artifact) || ~all(isfield(artifact, required)) || ...
            artifact.schemaVersion ~= 1 || ...
            ~isequaln(artifact.sweepIdentity, expectedIdentity) || ...
            ~sameExperiment(artifact.experimentSignature, expectedSignature)
        return;
    end
    payload = artifact.payload;
    reusable = true;
catch
    payload = [];
    reusable = false;
end
end

function path = normalizePath(value, label)
validType = (ischar(value) && isrow(value)) || ...
    (isstring(value) && isscalar(value));
if ~validType
    error('updateSweepCheckpoint:InvalidPath', ...
        '%s must be a character row or string scalar.', label);
end
path = string(value);
if ismissing(path) || strlength(path) == 0
    error('updateSweepCheckpoint:InvalidPath', ...
        '%s must not be empty.', label);
end
end

function name = artifactName(kind, target)
kind = string(validatestring(kind, ...
    {'linear','dense','pnnn','summary'}));
switch kind
    case "linear"
        name = 'linear_sweep.mat';
    case "dense"
        name = 'sweep_dense_source.mat';
    case "summary"
        name = 'complexity_sweep.mat';
    otherwise
        validateattributes(target, {'numeric'}, ...
            {'scalar','integer','positive','finite'});
        name = sprintf('pnnn_target_%04d.mat', target);
end
end

function valid = validExperimentSignature(value)
valid = isstruct(value) && isfield(value, 'digest') && ...
    strlength(string(value.digest)) > 0;
end

function value = sameExperiment(a, b)
value = validExperimentSignature(a) && validExperimentSignature(b) && ...
    string(a.digest) == string(b.digest) && ...
    isequaln(a.schemaVersion, b.schemaVersion);
end

function deleteIfPresent(filename)
filename = char(normalizePath(filename, 'temporary filename'));
if isfile(filename)
    delete(filename);
end
end
