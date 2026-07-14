function split = buildCommonComparisonSplit(x, y, cfg)
% buildCommonComparisonSplit - Build the shared full-signal protocol.
% The classic selector defines the 4% identification rows. An internal
% deterministic split is used only for hyperparameter selection.

x = x(:);
y = y(:);
if isempty(x) || numel(y) ~= numel(x) || ...
        any(~isfinite(x)) || any(~isfinite(y))
    error('buildCommonComparisonSplit:InvalidSignals', ...
        'x and y must be finite vectors of equal non-zero length.');
end
required = {'identificationFraction','identificationSeed', ...
    'internalTrainFraction','internalSplitSeed','amplitudeBinCount'};
if ~isstruct(cfg) || ~all(isfield(cfg, required))
    error('buildCommonComparisonSplit:InvalidConfig', ...
        'The split configuration is incomplete.');
end

previous_rng = rng;
cleanup = onCleanup(@() rng(previous_rng));
rng(cfg.identificationSeed, 'twister');
identification_indices = double(sel_indices( ...
    x, y, cfg.identificationFraction));
identification_indices = identification_indices(:);
expected_identification = floor(cfg.identificationFraction*numel(x));
if isempty(identification_indices) || ...
        numel(unique(identification_indices)) ~= ...
        numel(identification_indices) || ...
        any(identification_indices < 1) || ...
        any(identification_indices > numel(x)) || ...
        abs(numel(identification_indices) - expected_identification) > 2
    error('buildCommonComparisonSplit:InvalidIdentification', ...
        'The classic selector returned invalid identification rows.');
end
full_signal_indices = (1:numel(x)).';
assert(all(ismember(identification_indices, full_signal_indices)));
assert(numel(full_signal_indices) == numel(x));
assert(numel(unique(identification_indices)) == ...
    numel(identification_indices));
if isfield(cfg, 'expectedSignalLength') && ...
        isfield(cfg, 'expectedIdentificationSamples') && ...
        numel(x) == cfg.expectedSignalLength && ...
        numel(identification_indices) ~= cfg.expectedIdentificationSamples
    error('buildCommonComparisonSplit:UnexpectedIdentificationSize', ...
        'The configured capture must produce exactly %d identification rows.', ...
        cfg.expectedIdentificationSamples);
end

n_identification = numel(identification_indices);
n_bins = min(double(cfg.amplitudeBinCount), n_identification);
if n_bins < 1 || n_bins ~= floor(n_bins)
    error('buildCommonComparisonSplit:InvalidBinCount', ...
        'amplitudeBinCount must be a positive integer.');
end
train_fraction = double(cfg.internalTrainFraction);
if ~isscalar(train_fraction) || ~isfinite(train_fraction) || ...
        train_fraction <= 0 || train_fraction >= 1
    error('buildCommonComparisonSplit:InvalidTrainFraction', ...
        'internalTrainFraction must lie strictly between zero and one.');
end

[~, amplitude_order] = sort(abs(y(identification_indices)), 'ascend');
bin_edges = round(linspace(0, n_identification, n_bins + 1));
bin_counts = diff(bin_edges(:));
raw_train_counts = train_fraction*bin_counts;
train_counts = floor(raw_train_counts);
target_train_count = floor(train_fraction*n_identification);
remaining = target_train_count - sum(train_counts);
[~, remainder_order] = sort( ...
    raw_train_counts - train_counts, 'descend');
train_counts(remainder_order(1:remaining)) = ...
    train_counts(remainder_order(1:remaining)) + 1;

rng(cfg.internalSplitSeed, 'twister');
train_positions = false(n_identification, 1);
amplitude_bin = zeros(n_identification, 1);
for bin = 1:n_bins
    ordered_positions = amplitude_order( ...
        bin_edges(bin)+1:bin_edges(bin+1));
    amplitude_bin(ordered_positions) = bin;
    permutation = randperm(numel(ordered_positions));
    chosen = ordered_positions(permutation(1:train_counts(bin)));
    train_positions(chosen) = true;
end

internal_train_indices = sort(identification_indices(train_positions));
internal_validation_indices = sort(identification_indices(~train_positions));

if ~isempty(intersect(internal_train_indices, ...
        internal_validation_indices)) || ...
        ~isequal(sort([internal_train_indices; ...
        internal_validation_indices]), sort(identification_indices)) || ...
        numel(internal_train_indices) ~= target_train_count || ...
        ~all(ismember(identification_indices, full_signal_indices))
    error('buildCommonComparisonSplit:InvalidPartition', ...
        'The internal split does not partition the identification rows.');
end

split = struct();
split.internalTrainIndices = internal_train_indices;
split.internalValidationIndices = internal_validation_indices;
split.identificationIndices = identification_indices;
split.fullSignalIndices = full_signal_indices;
split.identificationAmplitudeBin = amplitude_bin;
split.identificationSeed = cfg.identificationSeed;
split.internalSplitSeed = cfg.internalSplitSeed;
split.identificationFraction = cfg.identificationFraction;
split.internalTrainFraction = train_fraction;
split.internalAmplitudeBinCount = n_bins;
split.signalLength = numel(x);
end
