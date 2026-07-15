function createComparisonPlots(results, resultDirectory)
% createComparisonPlots - Save the two established comparison figures.
% Figures use the frozen result table and do not perform model selection or
% scientific calculations beyond plotting the existing coordinates.

labels = string(results.Model);
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1150 680]);
scatter(results.NumRealParameters, results.FullSignalNMSEdB, 65, 'filled');
grid on;
xlabel('Active real parameters');
ylabel('Full-signal NMSE (dB)');
title('Full-signal NMSE versus active real parameters');
text(results.NumRealParameters, results.FullSignalNMSEdB, ...
    "  " + labels, 'Interpreter', 'none', 'FontSize', 8);
exportgraphics(figureHandle, fullfile(resultDirectory, ...
    'comparison_nmse_parameters.png'), 'Resolution', 160);
close(figureHandle);

markerSizes = 45 + 100*sqrt(results.NumRealParameters ./ ...
    max(results.NumRealParameters));
figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1150 680]);
scatter(results.FLOPsPerSample, results.FullSignalNMSEdB, ...
    markerSizes, 'filled');
grid on;
xlabel('FLOPs/sample');
ylabel('Full-signal NMSE (dB)');
title('Full-signal NMSE versus inference FLOPs/sample');
annotations = labels + " (p=" + string(results.NumRealParameters) + ")";
text(results.FLOPsPerSample, results.FullSignalNMSEdB, ...
    "  " + annotations, 'Interpreter', 'none', 'FontSize', 8);
exportgraphics(figureHandle, fullfile(resultDirectory, ...
    'comparison_flops.png'), 'Resolution', 160);
close(figureHandle);
end
