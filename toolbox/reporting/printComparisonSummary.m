function printComparisonSummary(study, split, resultDirectory)
% printComparisonSummary - Display the established final comparison tables.
% The summary reports frozen metrics only and states explicitly that the
% full-signal evaluation contains the identification rows.

disp(' ');
disp('=== Corrected full-signal main comparison ===');
disp(study.mainResults(:, {'Model','NumRealParameters', ...
    'IdentificationNMSEdB','FullSignalNMSEdB','FLOPsPerSample', ...
    'AdditionalOperationsPerSample'}));
disp(' ');
disp('=== All corrected-protocol models ===');
disp(study.comparisonResults(:, {'Model','NumRealParameters', ...
    'IdentificationNMSEdB','FullSignalNMSEdB','FLOPsPerSample'}));
fprintf('Complex/coupled full-signal relative error: %.6e\n', ...
    study.linear.equivalenceRelativeError);
fprintf('All final model fits used %d identification rows.\n', ...
    numel(split.identificationIndices));
fprintf(['Full-signal evaluation used %d rows and includes the ' ...
    'identification rows.\n'], numel(split.fullSignalIndices));
fprintf('Results: %s\n', resultDirectory);
end
