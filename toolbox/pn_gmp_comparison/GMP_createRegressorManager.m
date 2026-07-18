function rManagerGMP = GMP_createRegressorManager(x, y, cfg)
% GMP_createRegressorManager Builds the deterministic GMP basis without
% evaluating it on the signal. This avoids materializing the large GVG
% identification matrix in baseline-only workflows.

GVGconfig.Qpmax = cfg.Qpmax;
GVGconfig.Qnmax = cfg.Qnmax;
GVGconfig.Pmax  = cfg.Pmax;
GVGconfig.maxPopulation = cfg.maxPopulation;
GVGconfig.evaluationtype = 'maxPopulation';
GVGconfig.mutationrate = 0.7;
GVGconfig.crossoverrate = 0.5;
GVGconfig.verbosity = 0;
GVGconfig.showPlots = false;
GVGconfig.inittype = 'GMP';

rManagerGMP = regressorManager(x(:), y(:), GVGconfig);
rManagerGMP.initialization();
rManagerGMP.removerepeated();

nRegs = numel(rManagerGMP.regPopulation);
rManagerGMP.s = 1:nRegs;
rManagerGMP.nopt = nRegs;
rManagerGMP.nmse = NaN;
rManagerGMP.nmsev = [];
rManagerGMP.clearRegressors();
end
