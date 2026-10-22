function result = runPNIQAblation(sweepDirectory, P, savedPaths, makeFigures)
% Linear DPD ablation only: historical baselines, native IQ, and reciprocal support.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(root, 'toolbox')));
if nargin < 2, P = 310; end
if nargin < 3, savedPaths = {}; end % Global feature identifiers from the same capture/dictionary.
if nargin < 4, makeFigures = true; end
data = load(fullfile(sweepDirectory, 'experiment.mat'));
saved = load(fullfile(sweepDirectory, sprintf('models_P%04d.mat', P)), 'bundle');
cfg = data.cfg;
x = data.x; y = data.y; split = data.split; % Already DC-removed DPD u -> x_ILC.
models = saved.bundle.models(1:2);
historicalMetrics = saved.bundle.metrics(1:2,:);
clear saved;
cfg.sweep.parameterGrid = P;
manager = buildGMPManager(cfg.gmp);
population = (1:numel(manager.regPopulation)).';
K = numel(population);
descriptors = arrayfun(@(k) factorizeGMPRegressor(manager.regPopulation(k),k), population);
pnCandidates = find(~[false(K,1); [descriptors.QColumnStructurallyZero].']);
% Without PN only the imaginary envelope-only auxiliary is structurally zero.
iqCandidates = find(~[false(K,1); [descriptors.auxiliaryType].' == "envelope"]);
commonCandidates = intersect(iqCandidates, pnCandidates, 'stable');
labels = ["GMP historical", "PN-IQ historical", "IQ native DOMP", ...
    "IQ common DOMP", "PN-IQ on IQ common support"];
directory = fullfile(root, 'results', 'pniq_ablation', ...
    sprintf('P%d_%s', P, char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))));
assert(~isfolder(directory), 'Result directory already exists.');
mkdir(directory);
fprintf('Ablation output: %s\n',directory);
timer = tic;

%% Two real DOMP selections, sharing one identification matrix.
id = split.identificationIndices;
U = buildGMPRegressorRows(x,id,manager,population);
A = [real(U),imag(U)];
clear U;
candidateSets = {iqCandidates,commonCandidates};
for k = 1:2
    candidates = candidateSets{k};
    fprintf('%s: %d candidates, %d selected features...\n',labels(k+2),numel(candidates),P/2);
    if isempty(savedPaths)
        path = selectDOMPSupport(A(:,candidates),y(id),P/2,cfg.gmp.dompOptions.columnTolerance);
        support = candidates(path);
    else
        support = savedPaths{k}(1:P/2);
        assert(all(ismember(support,candidates)) && numel(unique(support)) == P/2);
    end
    assert(numel(support) == P/2);
    models{k+2} = fitIQ(A(:,support),y(id),support,manager,descriptors);
    fprintf('%s fitted (elapsed %.1f s).\n',labels(k+2),toc(timer));
end
clear A;

%% Reciprocal fixed-support control: all identifiers admissible in both spaces.
[present,path] = ismember(models{4}.support,pnCandidates);
assert(all(present));
pn = fitPNIQGMP(x,y,split,cfg,manager,population,path);
models{5} = pn.parameters;
clear pn;
rotation = ones(numel(id),1);
nonzero = abs(x(id)) ~= 0;
rotation(nonzero) = conj(x(id(nonzero)))./abs(x(id(nonzero)));
sources = unique(models{5}.featureMap.SourceRegressorIndex,'stable');
A = buildPNIQFeatures(x,id,rotation,manager,sources,descriptors);
A = A(:,models{5}.selectedColumns);
A = A./max(abs(A),[],1);
models{5}.normalizedRank = rank(A,max(size(A))*eps(norm(A,2)));
assert(models{5}.normalizedRank == P/2, 'Reciprocal PN support is rank deficient.');
clear A;
historicalPN = models{2}.featureMap;
historicalPNIds = historicalPN.SourceRegressorIndex + K*historicalPN.IsQ;
invalidTransfer = historicalPN(~ismember(historicalPNIds,iqCandidates),:);

%% Predictions and the same three NMSE domains; historical models are not fitted.
predictions = complex(zeros(numel(x),numel(models)));
auxiliary = cell(size(models));
metrics = table(); operations = cell(size(models));
for k = 1:numel(models)
    [predictions(:,k),auxiliary{k}] = predictPNIQAblation(models{k},x,cfg.gmp.blockSize);
    m = models{k};
    if m.kind == "complex"
        counts = countModelOperations(manager.regPopulation,m.support);
        counts = counts(1,:);
        parameters = 2*numel(m.coefficients);
    else
        counts = countModelOperations(manager.regPopulation, ...
            unique(m.featureMap.SourceRegressorIndex,'stable'),P/2);
        if m.kind == "iq"
            counts = countIQOperations(counts(1,:),m,descriptors);
        else
            counts = counts(4,:);
        end
        parameters = numel(m.coefficientsI)+numel(m.coefficientsQ);
    end
    assert(parameters == P);
    cost = countModelFLOPs(counts);
    operations{k} = counts;
    metrics = [metrics; table(labels(k),parameters, ...
        nmseComplexDb(y(id),predictions(id,k)), ...
        nmseComplexDb(y(split.evaluationIndices),predictions(split.evaluationIndices,k)), ...
        nmseComplexDb(y,predictions(:,k)),cost.FLOPsPerSample, ...
        'VariableNames',{'Model','RealParameters','IdentificationNMSEdB', ...
        'EvaluationNMSEdB','FullSignalNMSEdB','FLOPsPerSample'})]; %#ok<AGROW>
end
for field = ["IdentificationNMSEdB","EvaluationNMSEdB","FullSignalNMSEdB"]
    assert(max(abs(metrics.(field)(1:2)-historicalMetrics.(field))) < 1e-8);
    metrics.(field)(1:2) = historicalMetrics.(field); % Preserve the recorded baseline exactly.
end
assert(isequal(metrics.FLOPsPerSample(1:2),historicalMetrics.FLOPsPerSample));
result = struct('directory',directory,'sourceDirectory',sweepDirectory, ...
    'P',P,'config',cfg,'split',split,'models',{models},'metrics',metrics, ...
    'predictions',predictions,'operations',{operations}, ...
    'candidateCounts',[K,numel(iqCandidates),numel(pnCandidates),numel(commonCandidates)], ...
    'invalidHistoricalPNTransfer',invalidTransfer,'fitSeconds',toc(timer));
save(fullfile(directory,'ablation.mat'),'result','-v7.3');
writetable(metrics,fullfile(directory,'metrics.csv'));
disp(metrics);

%% Frozen-coefficient phase diagnostic: remove only the two auxiliary contributions.
capture = load(cfg.measurementFile,'x');
rawInput = double(capture.x(:));
rotations = table(); contributions = table();
types = ["conjugate","envelope"];
for k = 1:numel(models)
    reference = predictions(:,k);
    canonicalReference = reference - sum(auxiliary{k},2);
    for angle = [0 45 90 180 270]
        phase = exp(1j*angle*pi/180);
        if angle == 0
            prediction = reference; aux = auxiliary{k};
        else
            % Match ONLINE: rotate the original capture, then remove its DC.
            input = phase*rawInput;
            if cfg.pnnn.removeDC, input = input - mean(input); end
            [prediction,aux] = predictPNIQAblation(models{k},input,cfg.gmp.blockSize);
        end
        error = prediction-phase*reference;
        residual = (prediction-sum(aux,2))-phase*canonicalReference;
        rotations = [rotations; table(labels(k),angle,norm(error)/norm(reference), ...
            norm(residual)/norm(canonicalReference), ...
            norm(error-(sum(aux,2)-phase*sum(auxiliary{k},2)))/norm(reference), ...
            'VariableNames',{'Model','Degrees','EquivarianceError', ...
            'WithoutAuxiliariesError','UnexplainedErrorRelativeToBaseline'})]; %#ok<AGROW>
    end
    for a = 1:2
        contributions = [contributions; table(labels(k),types(a), ...
            norm(auxiliary{k}(:,a))/norm(reference), ...
            'VariableNames',{'Model','Auxiliary','RelativeOutputRMS'})]; %#ok<AGROW>
    end
end
result.rotations = rotations;
result.auxiliaryContributions = contributions;
result.totalSeconds = toc(timer);
save(fullfile(directory,'ablation.mat'),'result','-v7.3');
writetable(rotations,fullfile(directory,'rotations.csv'));
writetable(contributions,fullfile(directory,'auxiliary_contributions.csv'));
if makeFigures
s = getIEEEPaperStyle;
h = figure('Visible','off','Position',[100 100 850 470]);
plot(1:numel(labels),metrics.EvaluationNMSEdB,'o','MarkerFaceColor',s.gmpBlue,'MarkerSize',7);
xticks(1:numel(labels)); xticklabels(labels); xlim([0.5 numel(labels)+0.5]);
ylabel('Evaluation NMSE (dB)'); ylim([min(metrics.EvaluationNMSEdB)-0.5,max(metrics.EvaluationNMSEdB)+0.5]);
xtickangle(20); grid on; set(gca,'FontName',s.fontName,'FontSize',s.fontSize);
exportgraphics(h,fullfile(directory,'nmse_comparison.png'),'Resolution',200);
close(h);
at90 = rotations(rotations.Degrees == 90,:);
h = figure('Visible','off','Position',[100 100 750 440]);
semilogy((1:3)-0.04,at90.EquivarianceError(1:3),'o','Color',s.gmpBlue,'MarkerSize',7);
hold on;
semilogy((1:3)+0.04,at90.WithoutAuxiliariesError(1:3),'s','Color',s.pnOrange,'MarkerSize',7);
xticks(1:3); xticklabels(labels(1:3)); xlim([0.5 3.5]);
ylabel('Relative equivariance error at 90 deg'); grid on;
legend('Full predictor','Frozen auxiliary removal','Location','best');
set(gca,'FontName',s.fontName,'FontSize',s.fontSize);
exportgraphics(h,fullfile(directory,'equivariance_auxiliaries.png'),'Resolution',200);
close(h);
end
fprintf('Ablation complete: %.1f s.\n',result.totalSeconds);
end

function model = fitIQ(A,y,support,manager,descriptors)
% Same peak-scaled, unregularized lsqminnorm as the historical PN-IQ fit.
peaks = max(abs(A),[],1).';
assert(all(peaks > 0), 'Fixed support contains a zero column.');
B = A ./ peaks.'; outputPeak = max(abs(y));
tolerance = max(size(B))*eps(norm(B,2));
cI = lsqminnorm(B,real(y/outputPeak),tolerance)*outputPeak ./ peaks;
cQ = lsqminnorm(B,imag(y/outputPeak),tolerance)*outputPeak ./ peaks;
K = numel(manager.regPopulation);
source = mod(support-1,K)+1; isQ = support > K;
sources = unique(source,'stable');
[~,columns] = ismember(source,sources);
columns(isQ) = columns(isQ)+numel(sources);
model = struct('kind',"iq",'support',support, ...
    'regPopulation',compactGMPRegressors(manager,sources), ...
    'featureMap',table(source,isQ,'VariableNames',{'SourceRegressorIndex','IsQ'}), ...
    'descriptors',descriptors(sources),'selectedColumns',columns, ...
    'coefficientsI',cI,'coefficientsQ',cQ, ...
    'rankTolerance',tolerance,'normalizedRank',rank(B,tolerance));
assert(model.normalizedRank == numel(support), 'Selected IQ support is rank deficient.');
end

function counts = countIQOperations(counts,model,descriptors)
% Reuse magnitude/power costs, but form only the real components actually selected.
selected = descriptors(model.featureMap.SourceRegressorIndex);
nonlinear = [selected.canonicalGMP] & arrayfun(@(d) ~isempty(d.envelopeLags),selected).';
sources = unique(model.featureMap.SourceRegressorIndex,'stable');
uniqueDescriptors = descriptors(sources);
nNonlinear = nnz([uniqueDescriptors.canonicalGMP] & arrayfun(@(d) ~isempty(d.envelopeLags),uniqueDescriptors).');
F = height(model.featureMap);
counts.Model = "IQ-GMP";
counts.EffectiveRealFeatures = F; counts.RealParameters = 2*F;
counts.CoefficientMemoryBytes = 8*2*F;
counts.FeatureGenerationMultiplications = counts.FeatureGenerationMultiplications-2*nNonlinear+nnz(nonlinear);
counts.CoefficientMultiplications = 2*F; counts.CoefficientProductAdditions = 0;
counts.AccumulationAdditions = 2*(F-1);
counts.TotalRealMultiplicationsPerSample = counts.FeatureGenerationMultiplications+2*F;
counts.TotalRealAdditionsPerSample = counts.FeatureGenerationAdditions+2*(F-1);
end
