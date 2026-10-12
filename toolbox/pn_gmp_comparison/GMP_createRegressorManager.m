function rManagerGMP = GMP_createRegressorManager(x, y, cfg)
% GMP_createRegressorManager Builds the deterministic GMP basis without
% evaluating it on the signal. This avoids materializing the large GVG
% identification matrix in baseline-only workflows.

rManagerGMP = regressorManager(x(:), y(:), cfg);
rManagerGMP.initialization();
rManagerGMP.removerepeated();
end
