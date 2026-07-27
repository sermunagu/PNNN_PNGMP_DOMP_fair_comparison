clear;
clc;

%% Load the same configuration used by the simulation
projectRoot = fileparts(mfilename('fullpath'));

addpath(fullfile(projectRoot, 'config'));

for folder = ["complexity", "domp", "metrics", ...
        "pn_gmp_comparison", "pnnn", "splits", "sweep"]

    addpath(fullfile(projectRoot, 'toolbox', folder));
end

cfg = getFairDOMPComparisonConfig(projectRoot);

%% Load and prepare the same measurement
measurement = load(cfg.measurementFile, 'x', 'y');

[x, y] = selectXYByMapping( ...
    measurement.x, measurement.y, cfg.mappingMode);

x = x(:);
y = y(:);

if cfg.pnnn.removeDC
    x = x - mean(x);
    y = y - mean(y);
end

%% Reconstruct the common GMP candidate population
manager = GMP_createRegressorManager(x, y, cfg.gmp);

numComplexRegressors = numel(manager.regPopulation);

%% Count structurally zero PN-IQ Q-features
numStructuralZeroQFeatures = 0;

for index = 1:numComplexRegressors

    descriptor = factorizeGMPRegressor( ...
        manager.regPopulation(index), index);

    numStructuralZeroQFeatures = ...
        numStructuralZeroQFeatures + ...
        double(descriptor.QColumnStructurallyZero);
end

%% Dense parameter counts

% Complex GMP:
% one complex coefficient per complex regressor
complexGMPCoefficients = numComplexRegressors;
complexGMPRealParameters = 2 * numComplexRegressors;

% PN-IQ-GMP:
% each complex GMP regressor initially produces one I and one Q feature,
% except for structurally zero Q-features.
pniqCandidateFeatures = ...
    2 * numComplexRegressors - numStructuralZeroQFeatures;

% Each selected real feature has one coefficient for output I
% and another coefficient for output Q.
pniqRealParameters = 2 * pniqCandidateFeatures;

%% Display results
fprintf('\n=== Dense model sizes before sparse selection ===\n\n');

fprintf('Complex GMP:\n');
fprintf('  Complex candidate regressors: %d\n', ...
    complexGMPCoefficients);
fprintf('  Independent real parameters: %d\n\n', ...
    complexGMPRealParameters);

fprintf('PN-IQ-GMP:\n');
fprintf('  Candidate I/Q features: %d\n', ...
    pniqCandidateFeatures);
fprintf('  Structurally zero Q-features removed: %d\n', ...
    numStructuralZeroQFeatures);
fprintf('  Independent real parameters: %d\n\n', ...
    pniqRealParameters);

fprintf('Selected operating point:\n');
fprintf('  Retained real parameters: 340\n');