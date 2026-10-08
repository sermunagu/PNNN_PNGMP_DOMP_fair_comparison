function specs = compactGMPRegressors(manager, support)
% Save definitions only, not handle objects or cached regressor matrices.
specs = repmat(struct('X', [], 'Xconj', [], 'Xenv', []), numel(support), 1);
for k = 1:numel(support)
    regressor = manager.regPopulation(support(k));
    specs(k) = struct('X', regressor.X, 'Xconj', regressor.Xconj, ...
        'Xenv', regressor.Xenv);
end
end
