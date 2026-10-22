clear;
clc;

%% Load the same configuration used by the simulation
projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(fullfile(projectRoot, 'config'));

addpath(fullfile(projectRoot, 'toolbox', 'pn_gmp_comparison'));

cfg = getComparisonConfig(projectRoot);

%% Reconstruct the common GMP candidate population
% Population size depends on the dictionary configuration, not the capture.
manager = buildGMPManager(cfg.gmp);

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
    numComplexRegressors);
fprintf('  Independent real parameters: %d\n\n', ...
    complexGMPRealParameters);

fprintf('PN-IQ-GMP:\n');
fprintf('  Candidate I/Q features: %d\n', ...
    pniqCandidateFeatures);
fprintf('  Structurally zero Q-features removed: %d\n', ...
    numStructuralZeroQFeatures);
fprintf('  Independent real parameters: %d\n\n', ...
    pniqRealParameters);
