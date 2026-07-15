function layers = buildFairPNNNLayers(inputDimension, hiddenNeurons)
% buildFairPNNNLayers - One-hidden-layer sigmoid PNNN.

layers = [
    featureInputLayer(inputDimension, Name="input")
    fullyConnectedLayer(hiddenNeurons, Name="fc1")
    sigmoidLayer(Name="sigmoid1")
    fullyConnectedLayer(2, Name="fcOut")
];
end
