function count = countPNNNParameters(inputDimension, hiddenNeurons)
% countPNNNParameters - Count trainable real weights and biases in PNNN.
% The supported architecture has one real hidden layer and two real output
% channels; all trainable scalars are included exactly once.

validateattributes(inputDimension, {'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(hiddenNeurons, {'numeric'}, ...
    {'scalar','integer','positive','finite'});

count = struct();
count.inputDimension = double(inputDimension);
count.hiddenNeurons = double(hiddenNeurons);
count.realWeights = inputDimension*hiddenNeurons + 2*hiddenNeurons;
count.realBiases = hiddenNeurons + 2;
count.realParameters = count.realWeights + count.realBiases;
end
