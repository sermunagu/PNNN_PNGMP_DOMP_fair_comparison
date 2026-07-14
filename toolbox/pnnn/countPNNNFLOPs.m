function flops = countPNNNFLOPs(modelName, inputDimension, ...
    hiddenNeurons, memoryDepth, orders)
% countPNNNFLOPs - Count PNNN inference arithmetic under one explicit schedule.
% Phase rotation, reusable tap magnitudes, requested power chains, dense layers,
% biases, ELU calls, and final phase restoration are included per output row.

parameter_count = countPNNNParameters(inputDimension, hiddenNeurons);
orders = double(orders(:).');
if memoryDepth ~= 13 || ~isequal(orders, [1 3 5 7])
    error('countPNNNFLOPs:UnsupportedFeatureSchedule', ...
        'The audited FLOP schedule is defined for M=13 and orders [1 3 5 7].');
end
tap_count = memoryDepth + 1;
expected_dimension = 2*tap_count + tap_count*numel(orders);
if inputDimension ~= expected_dimension
    error('countPNNNFLOPs:DimensionMismatch', ...
        'The full PNNN feature dimension does not match M and orders.');
end

% Reuse |x(n)| from r(n), rotate all taps, compute the other magnitudes,
% and form powers 1/3/5/7 with a^2, a^3, a^5, a^7 (four M per tap).
phase_normalization_multiplications = 2;
phase_normalization_additions = 1;
feature_multiplications = 4*tap_count + 2*(tap_count-1) + 4*tap_count;
feature_additions = 2*tap_count + (tap_count-1);
network_multiplications = parameter_count.realWeights;
network_additions = parameter_count.realWeights;
phase_restoration_multiplications = 4;
phase_restoration_additions = 2;

Model = string(modelName);
NumRealParameters = parameter_count.realParameters;
RealMultiplicationsPerSample = phase_normalization_multiplications + ...
    feature_multiplications + network_multiplications + ...
    phase_restoration_multiplications;
RealAdditionsPerSample = phase_normalization_additions + ...
    feature_additions + network_additions + phase_restoration_additions;
ComplexMultiplicationsPerSample = 0;
ComplexAdditionsPerSample = 0;
NumELUPerSample = hiddenNeurons;
NumExpWorstCasePerSample = hiddenNeurons;
NumSqrtPerSample = tap_count;
NumRealDivisionsPerSample = 2;
NumAbsPerSample = tap_count;
PhaseNormalizationIncluded = true;
PhaseRestorationIncluded = true;
ComplexityScope = "inference per output sample";
NumRealWeights = parameter_count.realWeights;
NumRealBiases = parameter_count.realBiases;
DenseLayerCoreFLOPs = 2*inputDimension*hiddenNeurons;
OutputLayerCoreFLOPs = 4*hiddenNeurons;

specification = table(Model, NumRealParameters, ...
    RealMultiplicationsPerSample, RealAdditionsPerSample, ...
    ComplexMultiplicationsPerSample, ComplexAdditionsPerSample, ...
    NumELUPerSample, NumExpWorstCasePerSample, NumSqrtPerSample, ...
    NumRealDivisionsPerSample, NumAbsPerSample, ...
    PhaseNormalizationIncluded, PhaseRestorationIncluded, ...
    ComplexityScope, NumRealWeights, NumRealBiases, ...
    DenseLayerCoreFLOPs, OutputLayerCoreFLOPs);
flops = countModelFLOPs(specification, getFLOPConvention());
flops.DenseExecutionCoreFLOPsPerSample = flops.CoreFLOPsPerSample;
flops.IdealSparseCoreFLOPsPerSample = flops.CoreFLOPsPerSample;
flops.DenseMatrixFLOPsPerSample = flops.FLOPsPerSample;
flops.SparseZeroWeightsSkipped = false;
flops.IdealSparseRealMultiplicationsPerSample = ...
    flops.RealMultiplicationsPerSample;
flops.IdealSparseRealAdditionsPerSample = ...
    flops.RealAdditionsPerSample;
flops.IdealSparseCostRequiresSparseKernel = false;
end
