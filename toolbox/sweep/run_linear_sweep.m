function [sweep, models] = run_linear_sweep(x, y, split, cfg, paths)
% run_linear_sweep - Coordinate the two linear model sweeps.
% Complex GMP-DOMP and PN-IQ-GMP share one GMP population; each model
% owns one identification DOMP path, its prefix fits, predictions, and costs.

x = x(:);
y = y(:);
manager = GMP_createRegressorManager(x, y, cfg.gmp);

%% Report the initial dense Complex GMP size
Q = numel(manager.regPopulation);

denseComplexRegressors = Q;
denseComplexRealParameters = 2 * Q;

fprintf(['[Linear] Complex GMP dense population: %d complex ' ...
    'regressors, equivalent to %d real parameters.\n'], ...
    denseComplexRegressors, denseComplexRealParameters);

population = (1:Q).';

%% Fit the two scientific model families
if nargin < 5
    complexModel = fit_complex_gmp_domp(x, y, split, cfg, manager, population);
    pniqModel = fit_pniq_gmp(x, y, split, cfg, manager, population);
else
    % Explicit recovery: reuse the saved supports, never rerun DOMP.
    complexModel = fit_complex_gmp_domp(x, y, split, cfg, manager, population, paths.complex);
    pniqModel = fit_pniq_gmp(x, y, split, cfg, manager, population, paths.pniq);
end
models = {complexModel.parameters, pniqModel.parameters};




%% Package only the data reused after fitting
sweep.complexTable = complexModel.table;
sweep.pniqTable = pniqModel.table;
sweep.paths = struct('complex', complexModel.path, 'pniq', pniqModel.path);
sweep.predictions = struct('complexFull', complexModel.fullPredictions, ...
    'pniqFull', pniqModel.fullPredictions);
sweep.pniqPathMap = pniqModel.pniqPathMap;
sweep.coefficientRangeDefinition = ...
    string(cfg.sweep.coefficientRangeDefinition);
sweep.linearIdentificationScope = ...
    string(cfg.sweep.linearIdentificationScope);
sweep.linearPrincipalLambda = double(cfg.sweep.linearPrincipalLambda);
sweep.linearLambdaSelection = string(cfg.sweep.linearLambdaSelection);
sweep.fixedRidgeSupportPolicy = ...
    string(cfg.sweep.fixedRidgeSupportPolicy);
end
