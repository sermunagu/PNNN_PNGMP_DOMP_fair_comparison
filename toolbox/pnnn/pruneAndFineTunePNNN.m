function result = pruneAndFineTunePNNN(denseFit, features, neuralTargets, options)
% pruneAndFineTunePNNN - Prune one dense PNNN and fine-tune a fixed mask.
% Mode = "select-epochs" returns the best validation epoch.
% Mode = "fixed-epochs" runs the requested duration on all supplied training rows.

mode = string(options.Mode);
if ~isscalar(mode) || ~ismember(mode, ["select-epochs", "fixed-epochs"])
    error('pruneAndFineTunePNNN:InvalidMode', ...
        'Mode must be "select-epochs" or "fixed-epochs".');
end

trainingRows = options.TrainingRows;
targetActiveParameters = double(options.TargetActiveParameters);
config = options.Config;
if mode == "select-epochs"
    validationRows = options.ValidationRows;
else
    validationRows = [];
    config.pruning.fineTuneEpochs = double(options.FineTuneEpochs);
end

%% 1. Prune globally from the supplied immutable dense network
masks = createMagnitudePruningMasks(denseFit.network, targetActiveParameters);

sparseNetwork = applyLearnableMasks(denseFit.network, masks);
checkPruningMaskIntegrity(sparseNetwork, masks);

%% 2. Normalize train and optional validation rows with dense statistics
normalization = denseFit.normalization;
featuresTraining = (features(trainingRows, :) - normalization.muX) ./ normalization.sigmaX;
targetsTraining = (neuralTargets(trainingRows, :) - normalization.muY) ./ normalization.sigmaY;

if mode == "select-epochs"
    featuresValidation = (features(validationRows, :) - normalization.muX) ./ normalization.sigmaX;
    targetsValidation = (neuralTargets(validationRows, :) - normalization.muY) ./ normalization.sigmaY;
else
    featuresValidation = [];
    targetsValidation = [];
end

%% 3. Fine-tune while freezing every pruned weight at exactly zero
rng(double(options.NNSeed) + double(config.pruning.fineTuneSeedOffset), 'twister');
[sparseNetwork, bestEpoch] = fineTunePrunedNetwork( ...
    sparseNetwork, featuresTraining, targetsTraining, ...
    featuresValidation, targetsValidation, config, masks);

checkPruningMaskIntegrity(sparseNetwork, masks);

%% 4. Verify the exact active-parameter target and protected biases
counts = summarizeTrainableParameters(sparseNetwork, masks);

if counts.activeWeightParams + counts.activeBiasParams ~= targetActiveParameters
    error('pruneAndFineTunePNNN:TargetMismatch', 'Pruning did not produce the requested active parameter count.');
end

if counts.activeBiasParams ~= counts.totalBiasParams
    error('pruneAndFineTunePNNN:BiasProtectionFailure', 'Every bias must remain active and protected.');
end

result = struct('network', sparseNetwork, 'masks', {masks}, 'counts', counts, 'bestEpoch', bestEpoch);
end








function masks = createMagnitudePruningMasks( ...
    network, targetActiveParameters)
learnables = network.Learnables;
masks = cell(height(learnables), 1);
weightRows = zeros(0, 1);
weightSizes = cell(0, 1);
magnitudes = cell(0, 1);
totalParameters = 0;
totalWeights = 0;

for row = 1:height(learnables)
    value = learnableToNumeric(learnables.Value{row});
    masks{row} = true(size(value));
    totalParameters = totalParameters + numel(value);
    if lower(string(learnables.Parameter(row))) == "weights"
        weightRows(end+1, 1) = row; %#ok<AGROW>
        weightSizes{end+1, 1} = size(value); %#ok<AGROW>
        magnitudes{end+1, 1} = abs(value(:)); %#ok<AGROW>
        totalWeights = totalWeights + numel(value);
    end
end

protectedParameters = totalParameters - totalWeights;
if targetActiveParameters < protectedParameters || ...
        targetActiveParameters > totalParameters
    error('pruneAndFineTunePNNN:InvalidTarget', ...
        'The active-parameter target is incompatible with protected biases.');
end

pruneCount = totalParameters - targetActiveParameters;
prune = false(totalWeights, 1);
if pruneCount > 0
    allMagnitudes = vertcat(magnitudes{:});
    [~, order] = sort(allMagnitudes, 'ascend');
    prune(order(1:pruneCount)) = true;
end

offset = 0;
for index = 1:numel(weightRows)
    count = numel(magnitudes{index});
    localRows = offset + (1:count);
    masks{weightRows(index)} = reshape( ...
        ~prune(localRows), weightSizes{index});
    offset = offset + count;
end
end

function [bestNetwork, bestEpoch] = fineTunePrunedNetwork( ...
    network, featuresTraining, targetsTraining, ...
    featuresValidation, targetsValidation, config, masks)
trailingAverage = [];
trailingAverageSquared = [];
iteration = 0;
trainingCount = size(featuresTraining, 1);
useValidation = ~isempty(featuresValidation);
miniBatchSize = config.training.miniBatchSize;
dropPeriod = config.pruning.fineTuneLearnRateDropPeriod;
bestNetwork = network;
bestEpoch = 0;
bestValidationLoss = Inf;

for epoch = 1:config.pruning.fineTuneEpochs
    order = randperm(trainingCount);
    if config.training.verbose
        epochLoss = 0;
    end

    for first = 1:miniBatchSize:trainingCount
        batch = order(first:min(first + miniBatchSize - 1, trainingCount));
        dlX = dlarray(featuresTraining(batch, :).', "CB");
        dlT = dlarray(targetsTraining(batch, :).', "CB");
        iteration = iteration + 1;
        [loss, gradients] = dlfeval(@modelLoss, network, dlX, dlT);
        gradients = applyGradientMasks(gradients, masks);
        learnRate = config.pruning.fineTuneInitialLearnRate * ...
            config.training.learnRateDropFactor ^ ...
            floor((epoch-1) / max(dropPeriod, 1));
        [network, trailingAverage, trailingAverageSquared] = adamupdate( ...
            network, gradients, trailingAverage, trailingAverageSquared, ...
            iteration, learnRate);
        network = applyLearnableMasks(network, masks);
        if config.training.verbose
            epochLoss = epochLoss + ...
                double(gather(extractdata(loss))) * numel(batch);
        end
    end

    if useValidation
        validationLoss = computeValidationLoss(network, ...
            featuresValidation, targetsValidation, miniBatchSize);
        if isfinite(validationLoss) && validationLoss < bestValidationLoss
            bestValidationLoss = validationLoss;
            bestEpoch = epoch;
            bestNetwork = network;
        end
    else
        validationLoss = NaN;
        bestEpoch = epoch;
        bestNetwork = network;
    end

    if config.training.verbose
        trainLoss = epochLoss / trainingCount;
        fprintf(['Pruning fine-tune epoch %d/%d | train loss %.4g ' ...
            '| val loss %.4g\n'], epoch, ...
            config.pruning.fineTuneEpochs, trainLoss, validationLoss);
    end
end

if useValidation && bestEpoch == 0 && ...
        config.pruning.fineTuneEpochs > 0
    bestEpoch = config.pruning.fineTuneEpochs;
    bestNetwork = network;
end
bestNetwork = applyLearnableMasks(bestNetwork, masks);
end

function [loss, gradients] = modelLoss(network, dlX, dlT)
dlY = forward(network, dlX);
loss = mean((dlY - dlT).^2, "all");
gradients = dlgradient(loss, network.Learnables);
end

function value = computeValidationLoss( ...
    network, features, targets, miniBatchSize)
lossSum = 0;
count = 0;
for first = 1:miniBatchSize:size(features, 1)
    batch = first:min(first + miniBatchSize - 1, size(features, 1));
    dlX = dlarray(features(batch, :).', "CB");
    dlT = dlarray(targets(batch, :).', "CB");
    dlY = forward(network, dlX);
    loss = mean((dlY - dlT).^2, "all");
    lossSum = lossSum + ...
        double(gather(extractdata(loss))) * numel(batch);
    count = count + numel(batch);
end
value = lossSum / count;
end

function gradients = applyGradientMasks(gradients, masks)
for row = 1:height(gradients)
    gradients.Value{row} = applyMask(gradients.Value{row}, masks{row});
end
end

function network = applyLearnableMasks(network, masks)
if ~isa(network, 'dlnetwork')
    error('pruneAndFineTunePNNN:InvalidNetwork', ...
        'The network must be a dlnetwork object.');
end
learnables = network.Learnables;
if numel(masks) ~= height(learnables)
    error('pruneAndFineTunePNNN:MaskCountMismatch', ...
        'The mask count must match the network learnables.');
end
for row = 1:height(learnables)
    if ~isempty(masks{row})
        learnables.Value{row} = applyMask( ...
            learnables.Value{row}, masks{row});
    end
end
network.Learnables = learnables;
end

function value = applyMask(value, mask)
mask = logical(mask);
if isa(value, 'dlarray')
    data = extractdata(value);
else
    data = value;
end
maskValue = double(mask);
if isa(data, 'gpuArray')
    maskValue = gpuArray(maskValue);
end
maskValue = maskValue .* ones(size(data), "like", data);
value = value .* maskValue;
end

function checkPruningMaskIntegrity(network, masks)
learnables = network.Learnables;
for row = 1:height(learnables)
    value = learnableToNumeric(learnables.Value{row});
    if any(value(~logical(masks{row})) ~= 0, 'all')
        error('pruneAndFineTunePNNN:MaskViolation', ...
            'A pruned PNNN parameter became nonzero.');
    end
end
end

function data = learnableToNumeric(value)
if isa(value, 'dlarray')
    value = extractdata(value);
end
data = gather(value);
end
