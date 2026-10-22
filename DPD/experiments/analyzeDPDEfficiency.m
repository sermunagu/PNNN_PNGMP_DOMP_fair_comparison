function result = analyzeDPDEfficiency(sweepDirectory)
% Read the existing DPD sweep; no fitting or neural inference.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
if nargin < 1
    sweepDirectory = fullfile(root, 'results', 'dpd_sweep', 'tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748');
end
saved = load(fullfile(sweepDirectory, 'dpd_complexity_sweep.mat'), 'sweep');
sweep = saved.sweep;
data = load(fullfile(sweepDirectory, 'experiment.mat'), 'x', 'y', 'split', 'cfg');
cfg = data.cfg;
names = [cfg.names.complexGMPDOMP, cfg.names.pniqGMP, cfg.names.pnnn];
t = sweep.results;
assert(height(t) == 147 && isequal(unique(t.TargetRealParameters).', 20:10:500));
directory = fullfile(root, 'results', 'dpd_efficiency', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
assert(~isfolder(directory), 'Output folder already exists.');
mkdir(directory);
copyfile([mfilename('fullpath') '.m'], directory);

%% Pareto frontiers
% Strict Pareto dominance, with ties retained; smaller NMSE dB is better.
t.FamilyPareto = false(height(t), 1);
t.GlobalPareto = false(height(t), 1);
for i = 1:height(t)
    dominates = t.FLOPsPerSample <= t.FLOPsPerSample(i) & ...
        t.EvaluationNMSEdB <= t.EvaluationNMSEdB(i) & ...
        (t.FLOPsPerSample < t.FLOPsPerSample(i) | t.EvaluationNMSEdB < t.EvaluationNMSEdB(i));
    t.GlobalPareto(i) = ~any(dominates);
    t.FamilyPareto(i) = ~any(dominates & t.Model == t.Model(i));
end

%% Comparison at equal FLOPs
% Every observed cost is a ceiling. NaN means this family is unavailable.
ceilings = unique(t.FLOPsPerSample);
equalFLOPs = table();
for F = ceilings.'
    for name = names
        eligible = t(t.Model == name & t.FLOPsPerSample <= F, :);
        P = NaN;
        cost = NaN;
        nmse = NaN;
        if ~isempty(eligible)
            eligible = sortrows(eligible, {'EvaluationNMSEdB', 'FLOPsPerSample', 'TargetRealParameters'});
            P = eligible.TargetRealParameters(1);
            cost = eligible.FLOPsPerSample(1);
            nmse = eligible.EvaluationNMSEdB(1);
        end
        equalFLOPs = [equalFLOPs; table(F, name, P, cost, nmse, ...
            'VariableNames', {'CeilingFLOPs', 'Model', 'P', 'ActualFLOPs', 'EvaluationNMSEdB'})]; %#ok<AGROW>
    end
end

%% Comparison at equal NMSE
equalQuality = table();
for threshold = [-35 -37 -38 -38.5 -39]
    for name = names
        eligible = t(t.Model == name & t.EvaluationNMSEdB <= threshold, :);
        P = NaN;
        cost = NaN;
        nmse = NaN;
        if ~isempty(eligible)
            eligible = sortrows(eligible, {'FLOPsPerSample', 'TargetRealParameters', 'EvaluationNMSEdB'});
            P = eligible.TargetRealParameters(1);
            cost = eligible.FLOPsPerSample(1);
            nmse = eligible.EvaluationNMSEdB(1);
        end
        equalQuality = [equalQuality; table(threshold, name, ~isempty(eligible), P, cost, nmse, ...
            'VariableNames', {'TargetNMSEdB', 'Model', 'Reached', 'P', 'FLOPs', 'EvaluationNMSEdB'})]; %#ok<AGROW>
    end
end

%% Error by input amplitude
% Preserve each family's stored precision; do not concatenate single and double.
saved = load(fullfile(sweepDirectory, 'linear_sweep.mat'), 'linear');
k = find(saved.linear.complexTable.TargetRealParameters == 310);
predictions = {saved.linear.predictions.complexFull(:, k), saved.linear.predictions.pniqFull(:, k)};
clear saved;
saved = load(fullfile(sweepDirectory, 'pnnn_target_0310.mat'), 'point');
predictions{3} = saved.point.fullSignalPrediction;
evaluation = data.split.evaluationIndices;
assert(isequal(data.split.fullSignalIndices(:), (1:numel(data.x)).'));
assert(isequal(evaluation(:), setdiff(data.split.fullSignalIndices, data.split.identificationIndices)));
amplitude = abs(data.x(evaluation));
inputPeak = max(abs(data.x));
edges = unique(quantile(amplitude, linspace(0, 1, 21)));
bin = discretize(amplitude, edges);
referenceEnergy = sum(abs(data.y(evaluation)) .^ 2);
bins = table();
for family = 1:3
    expected = t.EvaluationNMSEdB(t.Model == names(family) & t.TargetRealParameters == 310);
    assert(abs(nmseComplexDb(data.y(evaluation), predictions{family}(evaluation)) - expected) < 1e-8);
    errorPower = abs(double(data.y(evaluation)) - double(predictions{family}(evaluation))) .^ 2;
    totalErrorEnergy = sum(errorPower);
    for b = 1:numel(edges) - 1
        rows = bin == b;
        targetEnergy = sum(abs(data.y(evaluation(rows))) .^ 2);
        errorEnergy = sum(errorPower(rows));
        meanAmplitude = mean(amplitude(rows));
        bins = [bins; table(names(family), b, nnz(rows), edges(b), edges(b + 1), ...
            meanAmplitude, meanAmplitude / inputPeak, targetEnergy, ...
            targetEnergy / referenceEnergy, errorEnergy, errorEnergy / nnz(rows), ...
            errorEnergy / totalErrorEnergy, 10 * log10(errorEnergy / targetEnergy), ...
            targetEnergy < 1e-8 * referenceEnergy, ...
            'VariableNames', {'Model', 'Bin', 'Samples', 'LowerAmplitude', 'UpperAmplitude', ...
            'MeanAmplitude', 'MeanNormalizedAmplitude', 'TargetEnergy', 'TargetEnergyFraction', ...
            'ErrorEnergy', 'MeanErrorPower', 'ErrorEnergyFraction', 'BinNMSEdB', 'LowReferenceEnergy'})]; %#ok<AGROW>
    end
end

%% Operation counts for the selected budget
saved = load(fullfile(sweepDirectory, 'models_P0310.mat'), 'bundle');
linearOperations = table();
costs = table();
for family = 1:2
    m = saved.bundle.models{family};
    ops = countModelOperations(m.regPopulation, (1:numel(m.regPopulation)).', 155);
    row = 1;
    if family == 2
        row = 4;
    end
    ops = ops(row, :);
    linearOperations = [linearOperations; ops]; %#ok<AGROW>
    cost = countModelFLOPs(ops);
    costs = [costs; cost]; %#ok<AGROW>
    assert(cost.FLOPsPerSample == t.FLOPsPerSample(t.Model == names(family) & t.TargetRealParameters == 310));
end
m = saved.bundle.models{3};
c = summarizeTrainableParameters(m.network, m.masks);
neuralCost = countSparsePNNNFLOPs(cfg.names.pnnn, 84, 12, 13, ...
    [1 3 5 7], c.activeWeightParams, c.activeBiasParams);
assert(neuralCost.FLOPsPerSample == t.FLOPsPerSample(t.Model == names(3) & t.TargetRealParameters == 310));

%% Save tables and export figures
result = struct('sourceDirectory', sweepDirectory, 'capture', cfg.measurementFile, ...
    'config', cfg, 'split', data.split, 'pareto', t, 'equalFLOPs', equalFLOPs, ...
    'equalQuality', equalQuality, 'amplitudeBins', bins, 'linearOperations', linearOperations, ...
    'linearCosts', costs, 'neuralCost', neuralCost, 'inputPeak', inputPeak, ...
    'amplitudeEdges', edges, 'predictions', {predictions}, 'directory', directory);
save(fullfile(directory, 'analysis.mat'), 'result', '-v7.3');
writetable(t, fullfile(directory, 'pareto.csv'));
writetable(equalFLOPs, fullfile(directory, 'equal_flops.csv'));
writetable(equalQuality, fullfile(directory, 'equal_quality.csv'));
writetable(bins, fullfile(directory, 'amplitude_bins_P310.csv'));
writetable(linearOperations, fullfile(directory, 'linear_operations_P310.csv'));
writetable(neuralCost, fullfile(directory, 'pnnn_cost_P310.csv'));

style = getIEEEPaperStyle;
addpath(cfg.paper.matlab2tikzSource);
h = figure('Visible', 'off', 'Position', [100 100 1000 420]);
tiledlayout(1, 2);
for panel = 1:2
    ax = nexttile;
    hold(ax, 'on');
    for family = 1:3
        rows = t.Model == names(family);
        scatter(ax, t.FLOPsPerSample(rows), t.EvaluationNMSEdB(rows), 15, ...
            style.mainColors(family, :), 'x', 'HandleVisibility', 'off');
        keep = rows & t.FamilyPareto;
        if panel == 2
            keep = rows & t.GlobalPareto;
        end
        q = sortrows(t(keep, :), 'FLOPsPerSample');
        plot(ax, q.FLOPsPerSample, q.EvaluationNMSEdB, 'o-', ...
            'Color', style.mainColors(family, :), 'DisplayName', names(family));
    end
    xlabel(ax, 'Arithmetic FLOPs/sample');
    ylabel(ax, 'Evaluation NMSE (dB)');
    titles = ["Within-family frontier", "Global frontier"];
    title(ax, titles(panel));
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end
exportPaperFigure(h, fullfile(directory, 'pareto'));
close(h);
h = figure('Visible', 'off');
hold on;
for family = 1:3
    rows = equalFLOPs.Model == names(family);
    stairs(equalFLOPs.CeilingFLOPs(rows), equalFLOPs.EvaluationNMSEdB(rows), ...
        'Color', style.mainColors(family, :), 'DisplayName', names(family), 'LineWidth', 1.6);
end
xlabel('FLOPs ceiling (observed grid)');
ylabel('Best feasible evaluation NMSE (dB)');
grid on;
legend('Location', 'best');
exportPaperFigure(h, fullfile(directory, 'equal_flops'));
close(h);
h = figure('Visible', 'off', 'Position', [100 100 1000 420]);
tiledlayout(1, 2);
for panel = 1:2
    ax = nexttile;
    hold(ax, 'on');
    for family = 1:3
        rows = bins.Model == names(family);
        values = bins.BinNMSEdB(rows);
        if panel == 2
            values = 100 * bins.ErrorEnergyFraction(rows);
        end
        plot(ax, bins.MeanNormalizedAmplitude(rows), values, 'o-', ...
            'Color', style.mainColors(family, :), 'DisplayName', names(family));
    end
    xlabel(ax, 'Mean |u-mean(u)| / full-record peak');
    labels = ["Within-bin NMSE (dB)", "Share of total error energy (%)"];
    ylabel(ax, labels(panel));
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end
exportPaperFigure(h, fullfile(directory, 'amplitude_error_P310'));
close(h);
disp(t(t.TargetRealParameters == 310, {'Model', 'EvaluationNMSEdB', 'FLOPsPerSample'}));
disp(equalQuality);
fprintf('Efficiency results: %s\n', directory);
end
