function rManagerGMP = GMP_createRegressorManager(x, y, cfg)
% GMP_createRegressorManager Builds the deterministic GMP basis without
% evaluating it on the signal. This avoids materializing the large GVG
% identification matrix in baseline-only workflows.

if nargin < 3 || isempty(cfg)
    cfg = struct();
end

GVGconfig.Qpmax = getCfgField(cfg, 'Qpmax', 50);
GVGconfig.Qnmax = getCfgField(cfg, 'Qnmax', 50);
GVGconfig.Pmax  = getCfgField(cfg, 'Pmax', 13);
GVGconfig.maxPopulation = getCfgField(cfg, 'maxPopulation', 100);
GVGconfig.evaluationtype = 'maxPopulation';
GVGconfig.mutationrate = 0.7;
GVGconfig.crossoverrate = 0.5;
GVGconfig.verbosity = getCfgField(cfg, 'verbosity', 0);
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

function value = getCfgField(cfg, name, defaultValue)
if isfield(cfg, name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
