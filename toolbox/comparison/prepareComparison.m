function [data, split, resultDirectory] = prepareComparison(cfg)
% prepareComparison - Load and prepare the shared modeled-block signals.
% The function applies the configured X/Y mapping and DC removal, validates
% the external data boundary, and builds the deterministic common split.

measurement = load(cfg.measurementFile, 'x', 'y');
if ~all(isfield(measurement, {'x', 'y'}))
    error('prepareComparison:InvalidMeasurement', ...
        'The configured measurement must contain x and y.');
end

[x, y] = selectXYByMapping( ...
    measurement.x, measurement.y, cfg.mappingMode);
x = x(:);
y = y(:);
if isempty(x) || numel(x) ~= numel(y) || ...
        any(~isfinite(x)) || any(~isfinite(y))
    error('prepareComparison:InvalidSignals', ...
        'Modeled-block X and Y must be aligned finite vectors.');
end

if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end

split = buildCommonComparisonSplit(x, y, cfg);
if ~isequal(split.fullSignalIndices(:), (1:numel(x)).') || ...
        ~all(ismember(split.identificationIndices, ...
        split.fullSignalIndices)) || ...
        ~isempty(intersect(split.internalTrainIndices, ...
        split.internalValidationIndices)) || ...
        ~isequal(sort([split.internalTrainIndices; ...
        split.internalValidationIndices]), ...
        sort(split.identificationIndices))
    error('prepareComparison:InvalidSplit', ...
        'The common identification and full-signal split is inconsistent.');
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
resultDirectory = fullfile(cfg.resultsRoot, timestamp);
if ~isfolder(resultDirectory)
    mkdir(resultDirectory);
end

data = struct('x', x, 'y', y);
end
