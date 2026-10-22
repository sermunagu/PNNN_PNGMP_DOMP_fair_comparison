function result = runGMPPhaseStudy(directory)
% Frozen diagnostics and fresh restricted DOMP. Resume only this explicit folder.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
source = fullfile(root, 'results', 'dpd_sweep', ...
    'tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748');
if nargin < 1
    directory = fullfile(root, 'results', 'gmp_phase_study', ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    assert(~isfolder(directory));
    mkdir(directory);
    copyfile([mfilename('fullpath') '.m'], directory);
end
if isfile(fullfile(directory, 'study.mat'))
    saved = load(fullfile(directory, 'study.mat'), 'result');
    result = saved.result;
    return;
end

%% Shared data and restricted dictionary
data = load(fullfile(source, 'experiment.mat'), 'x', 'y', 'split', 'cfg');
cfg = data.cfg;
x = data.x;
y = data.y;
split = data.split;
capture = load(fullfile(root, cfg.measurementFile), 'x');
raw = double(capture.x(:));
assert(isequal(raw - mean(raw), x), 'Preprocessing differs from the historical run.');
assert(isequal(split.fullSignalIndices(:), (1:numel(x)).'));
assert(isequal(split.evaluationIndices(:), ...
    setdiff(split.fullSignalIndices, split.identificationIndices)));
manager = buildGMPManager(cfg.gmp);
original = compactGMPRegressors(manager, (1:numel(manager.regPopulation)).');
q = arrayfun(@(r) numel(r.X) - numel(r.Xconj), original);
families = classifyTerms(original);
candidateIDs = find(q == 1);
classification = table((1:numel(q)).', families, q, q == 1, ...
    'VariableNames', {'OriginalRegressorIndex', 'Family', 'PhaseDegree', 'Retained'});
assert(numel(q) == 673 && numel(candidateIDs) == 671);
% Renumber the restricted population locally: helpers expect contiguous IDs.
restricted = struct('regPopulation', original(candidateIDs));
population = (1:numel(candidateIDs)).';
descriptors = arrayfun(@(k) factorizeGMPRegressor( ...
    restricted.regPopulation(k), k), population);
kept = find(~[false(numel(population), 1); [descriptors.QColumnStructurallyZero].']);
assert(numel(kept) == 1281);
id = split.identificationIndices;
ev = split.evaluationIndices;
rotation = ones(numel(id), 1);
nonzero = abs(x(id)) ~= 0;
rotation(nonzero) = conj(x(id(nonzero))) ./ abs(x(id(nonzero)));
angles = [0 30 60 90 180];
if ~isfile(fullfile(directory, 'classification.csv'))
    writetable(classification, fullfile(directory, 'classification.csv'));
    save(fullfile(directory, 'protocol.mat'), 'cfg', 'split', 'source', 'candidateIDs', 'angles');
end
fprintf('Phase study: %s\n', directory);

%% Frozen historical models first. Reuse prior 0/90/180-degree diagnostics.
frozenFile = fullfile(directory, 'frozen.mat');
if isfile(frozenFile)
    saved = load(frozenFile);
    frozen = saved.frozen;
    composition = saved.composition;
else
    frozen = table();
    composition = table();
    folders = ["P100_20261012_232926_441", "P310_20261012_223405_009", "P500_20261012_233006_611"];
    budgets = [100 310 500];
    for j = 1:3
        P = budgets(j);
        saved = load(fullfile(root, 'results', 'pniq_ablation', folders(j), 'ablation.mat'), 'result');
        previous = saved.result;
        bundle = load(fullfile(source, sprintf('models_P%04d.mat', P)), 'bundle');
        assert(isequaln(previous.models(1:2), bundle.bundle.models(1:2)));
        for family = 1:2
            model = previous.models{family};
            label = previous.metrics.Model(family);
            baseline = previous.predictions(:, family);
            baselineNorm = norm(baseline);
            for angle = angles
                old = previous.rotations.Model == label & previous.rotations.Degrees == angle;
                reused = any(old);
                if reused
                    error = previous.rotations.EquivarianceError(old);
                    without = previous.rotations.WithoutAuxiliariesError(old);
                else
                    phase = exp(1j * angle * pi / 180);
                    rotated = phase * raw;
                    rotated = rotated - mean(rotated);
                    prediction = predictPNIQAblation(model, rotated, cfg.gmp.blockSize);
                    error = norm(prediction - phase * baseline) / baselineNorm;
                    without = NaN;
                end
                frozen = [frozen; table(P, string(model.kind), angle, error, without, reused, ...
                    'VariableNames', {'P', 'Model', 'Degrees', 'RelativeError', 'FrozenAuxRemovalError', 'Reused'})]; %#ok<AGROW>
            end
            selectedFamilies = classifyTerms(model.regPopulation);
            for type = unique(families).'
                selected = selectedFamilies == type;
                features = nnz(selected);
                if model.kind == "pniq"
                    featureSources = mod(model.selectedColumns - 1, numel(model.regPopulation)) + 1;
                    features = nnz(selected(featureSources));
                end
                contribution = NaN;
                if type == "conjugate" || type == "envelope"
                    row = previous.auxiliaryContributions.Model == label & ...
                        previous.auxiliaryContributions.Auxiliary == type;
                    contribution = previous.auxiliaryContributions.RelativeOutputRMS(row);
                end
                composition = [composition; table(P, string(model.kind), type, nnz(families == type), ...
                    nnz(selected), features, contribution, 'VariableNames', ...
                    {'P', 'Model', 'Family', 'Candidates', 'SelectedSources', 'SelectedFeatures', 'RelativeContribution'})]; %#ok<AGROW>
            end
        end
        fprintf('Frozen P%d complete.\n', P);
    end
    save(frozenFile, 'frozen', 'composition');
    writetable(frozen, fullfile(directory, 'frozen_rotations.csv'));
    writetable(composition, fullfile(directory, 'historical_composition.csv'));
end

%% Two new paths to 250 atoms, shared by P310 first and P100/P500 afterwards.
for family = 1:2
    filename = fullfile(directory, sprintf('restricted_path_%d.mat', family));
    if isfile(filename)
        continue;
    end
    if family == 1
        A = buildGMPRegressorRows(x, id, restricted, population);
        target = y(id);
    else
        A = zeros(numel(id), numel(kept));
        for first = 1:cfg.sweep.candidateBlockSize:numel(id)
            rows = first:min(first + cfg.sweep.candidateBlockSize - 1, numel(id));
            block = buildPNIQFeatures(x, id(rows), rotation(rows), ...
                restricted, population, descriptors);
            A(rows, :) = block(:, kept);
        end
        target = rotation .* y(id);
    end
    fprintf('New restricted DOMP family %d, 250 atoms...\n', family);
    path = selectDOMPSupport(A, target, 250, cfg.gmp.dompOptions.columnTolerance);
    assert(numel(path) == 250);
    save(filename, 'path', 'candidateIDs');
    clear A block;
end

%% Validate P310 before extending; no old DOMP path is used.
metrics = table();
equivariance = table();
supports = table();
for P = [310 100 500]
    filename = fullfile(directory, sprintf('restricted_P%04d.mat', P));
    if isfile(filename)
        saved = load(filename, 'point');
        point = saved.point;
    else
        point = struct('P', P, 'models', {{}}, 'predictions', {{}}, ...
            'metrics', table(), 'equivariance', table(), 'supports', table());
        localCfg = cfg;
        localCfg.sweep.parameterGrid = P;
        historical = load(fullfile(source, sprintf('models_P%04d.mat', P)), 'bundle');
        for family = 1:2
            saved = load(fullfile(directory, sprintf('restricted_path_%d.mat', family)), 'path');
            path = saved.path(1:P / 2);
            if family == 1
                fit = fitComplexGMP(x, y, split, localCfg, restricted, population, path);
            else
                fit = fitPNIQGMP(x, y, split, localCfg, restricted, population, path);
            end
            model = fit.parameters;
            prediction = fit.fullPredictions;
            assert(fit.table.ActualRealParameters == P);
            assert(all(arrayfun(@(r) numel(r.X) - numel(r.Xconj) == 1, model.regPopulation)));
            point.models{family} = model;
            point.predictions{family} = prediction;
            if family == 1
                originalIDs = candidateIDs(path);
                isQ = false(P / 2, 1);
                oldIDs = historical.bundle.models{family}.support;
                A = buildGMPRegressorRows(x, id, restricted, path);
            else
                originalIDs = candidateIDs(model.featureMap.SourceRegressorIndex);
                isQ = model.featureMap.IsQ;
                old = historical.bundle.models{family}.featureMap;
                oldIDs = old.SourceRegressorIndex + numel(q) * old.IsQ;
                sources = (1:numel(model.regPopulation)).';
                A = buildPNIQFeatures(x, id, rotation, struct('regPopulation', model.regPopulation), ...
                    sources, model.descriptors);
                A = A(:, model.selectedColumns);
            end
            ids = originalIDs + numel(q) * isQ;
            A = A ./ max(abs(A), [], 1);
            singular = svd(A, 'econ');
            tolerance = max(size(A)) * eps(singular(1));
            rankValue = nnz(singular > tolerance);
            assert(rankValue == P / 2);
            evalNMSE = nmseComplexDb(y(ev), prediction(ev));
            baseline = historical.bundle.metrics(family, :);
            point.metrics = [point.metrics; table(P, string(model.kind), ...
                baseline.IdentificationNMSEdB, fit.table.IdentificationNMSEdB, ...
                baseline.EvaluationNMSEdB, evalNMSE, baseline.FullSignalNMSEdB, fit.table.FullSignalNMSEdB, ...
                baseline.FLOPsPerSample, fit.table.FLOPsPerSample, ...
                numel(intersect(ids, oldIDs)), rankValue, singular(1) / singular(end), ...
                'VariableNames', {'P', 'Model', 'HistoricalIdentificationNMSE', 'RestrictedIdentificationNMSE', ...
                'HistoricalEvaluationNMSE', 'RestrictedEvaluationNMSE', 'HistoricalFullNMSE', 'RestrictedFullNMSE', ...
                'HistoricalFLOPs', 'RestrictedFLOPs', 'CommonSupport', 'Rank', 'Condition2'})];
            point.supports = [point.supports; table(repmat(P, P / 2, 1), ...
                repmat(string(model.kind), P / 2, 1), (1:P / 2).', originalIDs, isQ, ...
                'VariableNames', {'P', 'Model', 'Order', 'OriginalRegressorIndex', 'IsQ'})];
            clear A;
            % Use the model's own replay as reference, not the blocked multi-fit matrix.
            replay = predictPNIQAblation(model, x, cfg.gmp.blockSize);
            replayError = norm(replay - prediction) / norm(prediction);
            assert(replayError < 1e-8, 'Stored model replay differs from fit.');
            replayNorm = norm(replay);
            for angle = angles
                error = 0;
                if angle ~= 0
                    phase = exp(1j * angle * pi / 180);
                    rotated = phase * raw;
                    rotated = rotated - mean(rotated);
                    shifted = predictPNIQAblation(model, rotated, cfg.gmp.blockSize);
                    error = norm(shifted - phase * replay) / replayNorm;
                end
                point.equivariance = [point.equivariance; table(P, string(model.kind), angle, error, replayError, ...
                    'VariableNames', {'P', 'Model', 'Degrees', 'RelativeError', 'ReplayRelativeError'})];
            end
        end
        save(filename, 'point', '-v7.3');
        disp(point.metrics);
        fprintf('Restricted P%d validated and saved.\n', P);
    end
    metrics = [metrics; point.metrics]; %#ok<AGROW>
    equivariance = [equivariance; point.equivariance]; %#ok<AGROW>
    supports = [supports; point.supports]; %#ok<AGROW>
end

%% Save the consolidated study and figures
result = struct('sourceDirectory', source, 'config', cfg, 'split', split, ...
    'classification', classification, 'composition', composition, 'frozen', frozen, ...
    'metrics', metrics, 'equivariance', equivariance, 'supports', supports, ...
    'minimumInputMagnitude', min(abs(x)), 'zeroInputSamples', nnz(x == 0), ...
    'inputPeak', max(abs(x)), 'directory', directory);
save(fullfile(directory, 'study.mat'), 'result', '-v7.3');
writetable(metrics, fullfile(directory, 'metrics.csv'));
writetable(equivariance, fullfile(directory, 'restricted_rotations.csv'));
writetable(supports, fullfile(directory, 'supports.csv'));
style = getIEEEPaperStyle;
addpath(cfg.paper.matlab2tikzSource);
h = figure('Visible', 'off');
hold on;
labels = ["complex", "pniq"];
for family = 1:2
    rows = metrics.Model == labels(family);
    t = sortrows(metrics(rows, :), 'P');
    plot(t.P, t.HistoricalEvaluationNMSE, 'o--', 'Color', style.mainColors(family, :), ...
        'DisplayName', labels(family) + " historical");
    plot(t.P, t.RestrictedEvaluationNMSE, 's-', 'Color', style.mainColors(family, :), ...
        'DisplayName', labels(family) + " restricted");
end
xlabel('Real parameters');
ylabel('Evaluation NMSE (dB)');
grid on;
legend('Location', 'best');
exportPaperFigure(h, fullfile(directory, 'restricted_nmse'));
close(h);
h = figure('Visible', 'off');
hold on;
for family = 1:2
    rows = frozen.P == 310 & frozen.Model == labels(family) & frozen.Degrees > 0;
    semilogy(frozen.Degrees(rows), frozen.RelativeError(rows), 'o--', ...
        'Color', style.mainColors(family, :), 'DisplayName', labels(family) + " historical");
    rows = equivariance.P == 310 & equivariance.Model == labels(family) & equivariance.Degrees > 0;
    semilogy(equivariance.Degrees(rows), equivariance.RelativeError(rows), 's-', ...
        'Color', style.mainColors(family, :), 'DisplayName', labels(family) + " restricted");
end
set(gca, 'YScale', 'log');
xlabel('Global input rotation (deg)');
ylabel('Relative equivariance error');
grid on;
legend('Location', 'best');
exportPaperFigure(h, fullfile(directory, 'phase_equivariance_P310'));
close(h);
fprintf('Phase study complete: %s\n', directory);
end

function family = classifyTerms(regressors)
family = strings(numel(regressors), 1);
for k = 1:numel(regressors)
    r = regressors(k);
    q = numel(r.X) - numel(r.Xconj);
    if q == -1
        family(k) = "conjugate";
    elseif q == 0
        family(k) = "envelope";
    elseif isempty(r.Xenv) && isequal(r.X, 0)
        family(k) = "current x";
    elseif isempty(r.Xenv)
        family(k) = "delayed x";
    else
        family(k) = "nonlinear GMP";
    end
end
end
