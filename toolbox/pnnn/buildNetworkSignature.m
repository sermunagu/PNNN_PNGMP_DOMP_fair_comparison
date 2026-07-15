function signature = buildNetworkSignature(denseFit)
% buildNetworkSignature - Hash one fitted dense PNNN and its normalization.
% The digest identifies the immutable dense source shared by every sparse
% sweep target; it does not replace the experiment/configuration signature.

if ~isstruct(denseFit) || ~isfield(denseFit, 'network') || ...
        ~isfield(denseFit, 'normalization') || ...
        ~isa(denseFit.network, 'dlnetwork')
    error('buildNetworkSignature:InvalidDenseFit', ...
        'A fitted dlnetwork and its normalization are required.');
end
requiredStats = {'muX','sigmaX','muY','sigmaY'};
if ~all(isfield(denseFit.normalization, requiredStats))
    error('buildNetworkSignature:InvalidNormalization', ...
        'The dense normalization is incomplete.');
end

digest = javaMethod('getInstance', ...
    'java.security.MessageDigest', 'SHA-256');
updateText(digest, "PNNN dense source signature v1");
learnables = denseFit.network.Learnables;
for row = 1:height(learnables)
    updateText(digest, string(learnables.Layer(row)));
    updateText(digest, string(learnables.Parameter(row)));
    value = learnables.Value{row};
    if isa(value, 'dlarray')
        value = extractdata(value);
    end
    updateNumeric(digest, gather(value));
end
for index = 1:numel(requiredStats)
    name = requiredStats{index};
    updateText(digest, name);
    updateNumeric(digest, denseFit.normalization.(name));
end
signature = struct('schemaVersion', 1, 'algorithm', "SHA-256", ...
    'digest', finishDigest(digest));
end

function updateNumeric(digest, value)
if ~isnumeric(value) || any(~isfinite(value), 'all')
    error('buildNetworkSignature:InvalidValue', ...
        'Dense learnables and normalization values must be finite numeric arrays.');
end
updateText(digest, class(value));
updateBytes(digest, numericBytes(uint64(size(value))));
updateBytes(digest, uint8(~isreal(value)));
updateBytes(digest, numericBytes(real(value(:))));
if ~isreal(value)
    updateBytes(digest, numericBytes(imag(value(:))));
end
end

function updateText(digest, value)
updateBytes(digest, unicode2native(char(string(value)), 'UTF-8'));
end

function bytes = numericBytes(value)
[~, ~, endian] = computer;
if endian == 'B'
    value = swapbytes(value);
end
bytes = typecast(value(:), 'uint8');
end

function updateBytes(digest, bytes)
digest.update(typecast(uint8(bytes(:)), 'int8'));
end

function value = finishDigest(digest)
bytes = typecast(int8(digest.digest()), 'uint8');
value = string(lower(reshape(dec2hex(bytes, 2).', 1, [])));
end
