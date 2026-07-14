function layers = buildFairPNNNLayers(inputDimension, hiddenNeurons)
% buildFairPNNNLayers - Build the fixed one-hidden-layer ELU PNNN.
% Inputs and outputs are real phase-normalized channels; the caller controls
% normalization, sample indices, initialization seed, and training protocol.

validateattributes(inputDimension, {'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(hiddenNeurons, {'numeric'}, ...
    {'scalar','integer','positive','finite'});

layers = [
    featureInputLayer(inputDimension, Name="input")
    fullyConnectedLayer(hiddenNeurons, Name="fc1")
    eluLayer(Name="elu1")
    fullyConnectedLayer(2, Name="fcOut")
];
end
