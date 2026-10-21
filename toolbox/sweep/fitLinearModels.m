function [linear, models] = fitLinearModels(x, y, split, cfg, paths)
% Fit both linear families for one or more budgets, optionally on saved supports.
% Complex GMP-DOMP and PN-IQ-GMP share one GMP population; each model
% owns one identification DOMP path, its prefix fits, predictions, and costs.

x = x(:);
y = y(:);
manager = buildGMPManager(x, y, cfg.gmp);

%% Report the initial dense Complex GMP size
Q = numel(manager.regPopulation);

fprintf(['[Linear] Complex GMP dense population: %d complex ' ...
    'regressors, equivalent to %d real parameters.\n'], ...
    Q, 2*Q);

population = (1:Q).';

%% Fit the two scientific model families
if nargin < 5
    complexModel = fitComplexGMP(x, y, split, cfg, manager, population);
    pniqModel = fitPNIQGMP(x, y, split, cfg, manager, population);
else
    % Explicit recovery: reuse the saved supports, never rerun DOMP.
    complexModel = fitComplexGMP(x, y, split, cfg, manager, population, paths.complex);
    pniqModel = fitPNIQGMP(x, y, split, cfg, manager, population, paths.pniq);
end
models = {complexModel.parameters, pniqModel.parameters};

%% Package only the data reused after fitting
linear.complexTable = complexModel.table;
linear.pniqTable = pniqModel.table;
linear.paths = struct('complex', complexModel.path, 'pniq', pniqModel.path);
linear.predictions = struct('complexFull', complexModel.fullPredictions, ...
    'pniqFull', pniqModel.fullPredictions);
linear.pniqPathMap = pniqModel.pniqPathMap;
linear.coefficientRangeDefinition = ...
    string(cfg.sweep.coefficientRangeDefinition);
linear.linearIdentificationScope = ...
    string(cfg.sweep.linearIdentificationScope);
linear.linearPrincipalLambda = double(cfg.sweep.linearPrincipalLambda);
linear.linearLambdaSelection = string(cfg.sweep.linearLambdaSelection);
linear.fixedRidgeSupportPolicy = ...
    string(cfg.sweep.fixedRidgeSupportPolicy);
end
