function manager = buildGMPManager(x, y, cfg)
% buildGMPManager Builds the deterministic GMP basis without
% evaluating it on the signal. This avoids materializing the large GVG
% identification matrix in baseline-only workflows.

manager = regressorManager(x(:), y(:), cfg);
manager.buildPopulation();
manager.finalizePopulation();
end
