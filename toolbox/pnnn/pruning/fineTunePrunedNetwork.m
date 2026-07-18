function [bestNet, bestEpoch] = fineTunePrunedNetwork( ...
    net, inputTrain, targetTrain, inputValidation, targetValidation, ...
    cfg, masks)
% fineTunePrunedNetwork - Fine-tune a PNNN while keeping pruned weights zero.
% Validation selects the best epoch; an empty validation set returns the
% fixed-epoch final network used by the final identification refit.

trailingAverage = [];
trailingAverageSquared = [];
iteration = 0;
trainingCount = size(inputTrain, 1);
useValidation = ~isempty(inputValidation);
miniBatchSize = cfg.training.miniBatchSize;
dropPeriod = cfg.pruning.fineTuneLearnRateDropPeriod;
bestNet = net;
bestEpoch = 0;
bestValidationLoss = Inf;

for epoch = 1:cfg.pruning.fineTuneEpochs
    order = randperm(trainingCount);
    if cfg.training.verbose
        epochLoss = 0;
    end

    for first = 1:miniBatchSize:trainingCount
        batch = order(first:min(first + miniBatchSize - 1, trainingCount));
        dlX = dlarray(inputTrain(batch, :).', "CB");
        dlT = dlarray(targetTrain(batch, :).', "CB");
        iteration = iteration + 1;
        [loss, gradients] = dlfeval(@modelLoss, net, dlX, dlT);
        gradients = applyGradientMasks(gradients, masks);
        learnRate = cfg.pruning.fineTuneInitialLearnRate * ...
            cfg.training.learnRateDropFactor ^ ...
            floor((epoch-1) / max(dropPeriod, 1));
        [net, trailingAverage, trailingAverageSquared] = adamupdate( ...
            net, gradients, trailingAverage, trailingAverageSquared, ...
            iteration, learnRate);
        net = applyLearnableMasks(net, masks);
        if cfg.training.verbose
            epochLoss = epochLoss + ...
                double(gather(extractdata(loss))) * numel(batch);
        end
    end

    if useValidation
        validationLoss = computeLoss( ...
            net, inputValidation, targetValidation, miniBatchSize);
        if isfinite(validationLoss) && validationLoss < bestValidationLoss
            bestValidationLoss = validationLoss;
            bestEpoch = epoch;
            bestNet = net;
        end
    else
        validationLoss = NaN;
        bestEpoch = epoch;
        bestNet = net;
    end

    if cfg.training.verbose
        trainLoss = epochLoss / trainingCount;
        fprintf("Pruning fine-tune epoch %d/%d | train loss %.4g | val loss %.4g\n", ...
            epoch, cfg.pruning.fineTuneEpochs, trainLoss, validationLoss);
    end
end

if useValidation && bestEpoch == 0 && cfg.pruning.fineTuneEpochs > 0
    bestEpoch = cfg.pruning.fineTuneEpochs;
    bestNet = net;
end
bestNet = applyLearnableMasks(bestNet, masks);
end

function [loss, gradients] = modelLoss(net, dlX, dlT)
dlY = forward(net, dlX);
loss = mean((dlY - dlT).^2, "all");
gradients = dlgradient(loss, net.Learnables);
end

function value = computeLoss(net, input, target, miniBatchSize)
lossSum = 0;
count = 0;
for first = 1:miniBatchSize:size(input, 1)
    batch = first:min(first + miniBatchSize - 1, size(input, 1));
    dlX = dlarray(input(batch, :).', "CB");
    dlT = dlarray(target(batch, :).', "CB");
    dlY = forward(net, dlX);
    loss = mean((dlY - dlT).^2, "all");
    lossSum = lossSum + double(gather(extractdata(loss))) * numel(batch);
    count = count + numel(batch);
end
value = lossSum / count;
end

function gradients = applyGradientMasks(gradients, masks)
for row = 1:height(gradients)
    gradients.Value{row} = applyMask(gradients.Value{row}, masks{row});
end
end

function value = applyMask(value, mask)
data = extractdata(value);
maskValue = double(mask) .* ones(size(data), "like", data);
if isa(data, 'gpuArray')
    maskValue = gpuArray(maskValue);
end
value = value .* maskValue;
end
