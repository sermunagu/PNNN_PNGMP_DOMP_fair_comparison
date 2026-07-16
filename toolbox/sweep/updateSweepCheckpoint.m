function [payload, reusable, filename] = updateSweepCheckpoint( ...
    resultDirectory, artifactKind, target, sweepIdentity, ...
    experimentSignature, newPayload)
% updateSweepCheckpoint - Load or atomically replace one sweep artifact.
% Artifact names are deterministic. Each call reads only the requested unit,
% so resuming one PNNN target never scans or rewrites previous targets.

if nargin < 6
    newPayload = [];
end
validateattributes(resultDirectory, {'char','string'}, {'scalartext'});
if ~isstruct(sweepIdentity) || ~isfield(sweepIdentity, 'digest') || ...
        ~validExperimentSignature(experimentSignature)
    error('updateSweepCheckpoint:InvalidIdentity', ...
        'Signed experiment and sweep identities are required.');
end
if ~isfolder(resultDirectory)
    mkdir(resultDirectory);
end
filename = fullfile(resultDirectory, artifactName(artifactKind, target));

if isempty(newPayload)
    [payload, reusable] = loadArtifact( ...
        filename, sweepIdentity, experimentSignature);
    return;
end

checkpointArtifact = struct('schemaVersion', 1, ...
    'sweepIdentity', sweepIdentity, ...
    'experimentSignature', experimentSignature, ...
    'payload', newPayload);
[directory, base, extension] = fileparts(filename);
temporary = fullfile(directory, [base '.tmp' extension]);
if isfile(temporary)
    delete(temporary);
end
cleanup = onCleanup(@() deleteIfPresent(temporary));
save(temporary, 'checkpointArtifact', '-v7.3');
[verified, reusable] = loadArtifact( ...
    temporary, sweepIdentity, experimentSignature);
if ~reusable || isempty(verified)
    error('updateSweepCheckpoint:SerializationFailed', ...
        'The temporary sweep artifact failed verification.');
end
movefile(temporary, filename, 'f');
payload = newPayload;
reusable = true;
clear cleanup;
end

function [payload, reusable] = loadArtifact( ...
    filename, expectedIdentity, expectedSignature)
payload = [];
reusable = false;
if ~isfile(filename)
    return;
end
try
    saved = load(filename, 'checkpointArtifact');
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
if isfile(filename)
    delete(filename);
end
end
