function split = buildCommonComparisonSplit(x, y, cfg)
% buildCommonComparisonSplit - Shared identification/train/validation rows.

x = x(:);
y = y(:);
n_samples = numel(x);

% Keep the original deterministic identification selector.
previous_rng = rng;
cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(cfg.identificationSeed, 'twister');
identification = double(sel_indices(x, y, cfg.identificationFraction));
identification = identification(:);
full_signal = (1:n_samples).';

% Split identification into amplitude-balanced train and validation subsets.
n_identification = numel(identification);
n_bins = min(cfg.amplitudeBinCount, n_identification);
[~, amplitude_order] = sort(abs(y(identification)), 'ascend');
bin_edges = round(linspace(0, n_identification, n_bins + 1));
bin_counts = diff(bin_edges(:));

raw_train_counts = cfg.internalTrainFraction*bin_counts;
train_counts = floor(raw_train_counts);
remaining = floor(cfg.internalTrainFraction*n_identification) - sum(train_counts);
[~, order] = sort(raw_train_counts - train_counts, 'descend');
train_counts(order(1:remaining)) = train_counts(order(1:remaining)) + 1;

rng(cfg.internalSplitSeed, 'twister');
is_train = false(n_identification, 1);
amplitude_bin = zeros(n_identification, 1);
for bin = 1:n_bins
    positions = amplitude_order(bin_edges(bin)+1:bin_edges(bin+1));
    amplitude_bin(positions) = bin;
    positions = positions(randperm(numel(positions)));
    is_train(positions(1:train_counts(bin))) = true;
end

split = struct( ...
    'internalTrainIndices', sort(identification(is_train)), ...
    'internalValidationIndices', sort(identification(~is_train)), ...
    'identificationIndices', identification, ...
    'fullSignalIndices', full_signal, ...
    'identificationAmplitudeBin', amplitude_bin, ...
    'identificationSeed', cfg.identificationSeed, ...
    'internalSplitSeed', cfg.internalSplitSeed, ...
    'identificationFraction', cfg.identificationFraction, ...
    'internalTrainFraction', cfg.internalTrainFraction, ...
    'internalAmplitudeBinCount', n_bins, ...
    'signalLength', n_samples);
end
