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
if activeBiases ~= dense_count.realBiases
    error('countSparsePNNNFLOPs:BiasCountMismatch', ...
        'The audited sparse schedule protects every bias.');
end

ideal_sparse_multiplications = dense_flops.RealMultiplicationsPerSample - ...
    dense_count.realWeights + activeWeights;
ideal_sparse_additions = dense_flops.RealAdditionsPerSample - ...
    dense_count.realWeights + activeWeights;

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
