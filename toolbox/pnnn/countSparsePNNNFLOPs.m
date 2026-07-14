function flops = countSparsePNNNFLOPs(modelName, inputDimension, ...
    hiddenNeurons, memoryDepth, orders, activeWeights, activeBiases)
% countSparsePNNNFLOPs - Separate dense and ideal sparse PNNN core costs.
% The primary FLOPs/sample count skips exactly zero weights. A secondary
% dense-matrix count is retained for implementations that do not skip zeros.

dense_flops = countPNNNFLOPs(modelName, inputDimension, ...
    hiddenNeurons, memoryDepth, orders);
dense_count = countPNNNParameters(inputDimension, hiddenNeurons);
validateattributes(activeWeights, {'numeric'}, ...
    {'scalar','integer','nonnegative','finite','<=',dense_count.realWeights});
validateattributes(activeBiases, {'numeric'}, ...
    {'scalar','integer','nonnegative','finite','<=',dense_count.realBiases});
if activeBiases ~= dense_count.realBiases
    error('countSparsePNNNFLOPs:BiasCountMismatch', ...
        'The audited sparse schedule protects every bias.');
end

dense_network_multiplications = dense_count.realWeights;
dense_network_additions = dense_count.realWeights;
non_network_multiplications = ...
    dense_flops.RealMultiplicationsPerSample - ...
    dense_network_multiplications;
non_network_additions = dense_flops.RealAdditionsPerSample - ...
    dense_network_additions;
ideal_sparse_multiplications = non_network_multiplications + activeWeights;
ideal_sparse_additions = non_network_additions + activeWeights;

flops = dense_flops;
flops.NumRealParameters = activeWeights + activeBiases;
flops.NumRealWeights = activeWeights;
flops.NumRealBiases = activeBiases;
flops.DenseExecutionCoreFLOPsPerSample = ...
    dense_flops.CoreFLOPsPerSample;
flops.IdealSparseCoreFLOPsPerSample = ...
    ideal_sparse_multiplications + ideal_sparse_additions;
flops.IdealSparseRealMultiplicationsPerSample = ...
    ideal_sparse_multiplications;
flops.IdealSparseRealAdditionsPerSample = ideal_sparse_additions;
flops.IdealSparseCostRequiresSparseKernel = true;
flops.FLOPsPerSample = ideal_sparse_multiplications + ...
    ideal_sparse_additions;
flops.DenseMatrixFLOPsPerSample = ...
    dense_flops.FLOPsPerSample;
flops.SparseZeroWeightsSkipped = true;
end
