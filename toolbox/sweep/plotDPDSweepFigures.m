function files = plotDPDSweepFigures(sweep)
% Export saved ILC-imitation metrics only; no model loading or fitting.
cfg = sweep.config;
addpath(cfg.paper.matlab2tikzSource);
directory = char(sweep.resultDirectory);
% Short display name keeps the nine-entry Ridge legend inside the figure.
results = sweep.results;
fixed = sweep.fixedResults;
results.Model(results.Model == cfg.names.complexGMPDOMP) = "GMP-DOMP";
fixed.Model(fixed.Model == cfg.names.complexGMPDOMP) = "GMP-DOMP";
cfg.names.complexGMPDOMP = "GMP-DOMP";
options = struct('names', cfg.names, 'fixedLambdas', cfg.fixedRidgeLambdas, ...
    'isNMSE', false); % Autoscale: the historical PA -30 dB cap is inappropriate here.
variables = ["IdentificationNMSEdB", "EvaluationNMSEdB", "FullSignalNMSEdB"];
labels = ["ILC NMSE - identification (dB)", ...
    "ILC NMSE - evaluation complement (dB)", ...
    "ILC NMSE - full signal (dB)"];
keys = ["identification", "evaluation", "full_signal"];
files = struct();
for index = 1:numel(variables)
    options.metricVariable = char(variables(index));
    options.metricLabel = char(labels(index));
    options.includeFixed = true;
    files.(keys(index) + "Parameters") = plotSweepPaperFigure( ...
        results, fixed, 'ActualRealParameters', ...
        'Active real parameters', fullfile(directory, "dpd_nmse_" + keys(index) + "_parameters"), options);
    options.includeFixed = false;
    files.(keys(index) + "FLOPs") = plotSweepPaperFigure( ...
        results, table(), 'FLOPsPerSample', 'FLOPs per sample', ...
        fullfile(directory, "dpd_nmse_" + keys(index) + "_flops"), options);
end
options.metricVariable = 'FLOPsPerSample';
options.metricLabel = 'FLOPs per sample (PNNN: ideal sparse)';
files.flopsParameters = plotSweepPaperFigure(results, table(), ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(directory, 'dpd_flops_parameters'), options);
options.metricVariable = 'MaxAbsRealParameter';
options.metricLabel = 'Maximum absolute coefficient/parameter';
options.includeFixed = true;
options.useLogWhenPositive = true;
files.maxAbsParameter = plotSweepPaperFigure(results, fixed, ...
    'ActualRealParameters', 'Active real parameters', ...
    fullfile(directory, 'dpd_max_abs_parameter'), options);

style = getIEEEPaperStyle();
h = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', 'Position', [1 1 7 4.5]);
cleanup = onCleanup(@() close(h));
ax = axes(h);
rows = sweep.results.Model == cfg.names.pnnn;
plot(ax, sweep.results.ActualRealParameters(rows), sweep.results.WeightSparsityPercent(rows), ...
    'Color', style.pnnnGreen, 'LineStyle', style.mainLineStyles{3}, ...
    'Marker', style.mainMarkers{3}, 'MarkerSize', style.markerSize, 'LineWidth', style.mainLineWidth);
xlabel(ax, 'Active real parameters');
ylabel(ax, 'PNNN weight sparsity (%)');
set(ax, 'FontName', style.fontName, 'FontSize', style.fontSize, 'LineWidth', style.axisLineWidth);
grid(ax, 'on');
box(ax, 'on');
files.sparsityParameters = exportPaperFigure(h, fullfile(directory, 'dpd_pnnn_sparsity_parameters'));
end
