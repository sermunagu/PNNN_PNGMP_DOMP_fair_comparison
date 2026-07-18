function sweep = run_linear_sweep(x, y, split, cfg)
% run_linear_sweep - Coordinate the two linear model sweeps.
% Complex GMP and independent PN-IQ share one GMP population and split;
% each model owns its DOMP paths, fits, predictions, metrics, and costs.

x = x(:);
y = y(:);
manager = GMP_createRegressorManager(x, y, cfg.gmp);
population = (1:numel(manager.regPopulation)).';

%% Fit the two scientific model families
complexModel = fit_complex_gmp_domp(x, y, split, cfg, manager, population);
pnModel = fit_independent_pniq_domp(x, y, split, cfg, manager, population);




%% Package only the data reused after fitting
sweep.complexTable = complexModel.table;
sweep.pnTable = pnModel.table;
sweep.paths = struct('complex', complexModel.path, 'pn', pnModel.path);
sweep.predictions = struct( ...
    'complexFull', complexModel.fullPredictions, ...
    'pnFull', pnModel.fullPredictions);
sweep.pnPathMap = pnModel.pnPathMap;
end
