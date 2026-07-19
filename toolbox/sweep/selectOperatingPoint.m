function selection = selectOperatingPoint(results, selectionConfig)
% selectOperatingPoint - Quantify and select the sparse PNNN operating point.
% The admissible set is near-optimal PNNN points that are no worse than the
% Complex GMP point at the same real-parameter budget. Complexity breaks ties.

if nargin < 2 || isempty(selectionConfig)
    selectionConfig = struct();
end
if ~isfield(selectionConfig, 'nmseToleranceDb')
    selectionConfig.nmseToleranceDb = 0.20;
end
if ~isfield(selectionConfig, 'sensitivityTolerancesDb')
    selectionConfig.sensitivityTolerancesDb = [0.10 0.15 0.20 0.25];
end
criterionName = "near-optimal minimum-complexity criterion";
required = ["Model", "ActualRealParameters", "FullSignalNMSEdB", ...
    "FLOPsPerSample"];
if ~istable(results) || ~all(ismember(required, ...
        string(results.Properties.VariableNames)))
    error('selectOperatingPoint:InvalidResults', ...
        'The canonical results table is missing required variables.');
end

pnnnMask = string(results.Model) == "Sparse PNNN N12";
gmpMask = string(results.Model) == "Complex GMP DOMP sweep";
pnnn = results(pnnnMask, :);
gmp = results(gmpMask, :);
if isempty(pnnn) || isempty(gmp)
    error('selectOperatingPoint:MissingFamily', ...
        'Both sparse PNNN and Complex GMP sweep rows are required.');
end
[~, order] = sortrows([double(pnnn.ActualRealParameters), ...
    double(pnnn.FLOPsPerSample)], [1 2]);
pnnn = pnnn(order, :);
if numel(unique(pnnn.ActualRealParameters)) ~= height(pnnn)
    error('selectOperatingPoint:DuplicatePNNNBudget', ...
        'Each PNNN real-parameter budget must occur exactly once.');
end

gmpNMSE = zeros(height(pnnn), 1);
for row = 1:height(pnnn)
    match = gmp.ActualRealParameters == pnnn.ActualRealParameters(row);
    if nnz(match) ~= 1
        error('selectOperatingPoint:MissingGMPPair', ...
            'PNNN budget %g does not have exactly one GMP pair.', ...
            pnnn.ActualRealParameters(row));
    end
    gmpNMSE(row) = gmp.FullSignalNMSEdB(match);
end

pnnnNMSE = double(pnnn.FullSignalNMSEdB);
pnnnFLOPs = double(pnnn.FLOPsPerSample);
pnnnParameters = double(pnnn.ActualRealParameters);
if any(~isfinite([pnnnNMSE; pnnnFLOPs; pnnnParameters; gmpNMSE]))
    error('selectOperatingPoint:NonfiniteData', ...
        'Selection metrics must be finite.');
end
bestNMSE = min(pnnnNMSE);
bestCandidateRows = find(pnnnNMSE == bestNMSE);
[~, localBest] = sortrows([pnnnFLOPs(bestCandidateRows), ...
    pnnnParameters(bestCandidateRows)], [1 2]);
bestRow = bestCandidateRows(localBest(1));

NMSELossFromBestDb = pnnnNMSE - bestNMSE;
PNNNGainOverGMPDb = gmpNMSE - pnnnNMSE;
FLOPsSavedVsBest = pnnnFLOPs(bestRow) - pnnnFLOPs;
FLOPsSavedPercentVsBest = 100*FLOPsSavedVsBest/pnnnFLOPs(bestRow);
ParametersSavedVsBest = pnnnParameters(bestRow) - pnnnParameters;
ParametersSavedPercentVsBest = ...
    100*ParametersSavedVsBest/pnnnParameters(bestRow);
GainFromPreviousPointDb = nan(height(pnnn), 1);
ExtraFLOPsFromPreviousPoint = nan(height(pnnn), 1);
MarginalGainDbPer100AdditionalFLOPs = nan(height(pnnn), 1);
if height(pnnn) > 1
    GainFromPreviousPointDb(2:end) = pnnnNMSE(1:end-1) - pnnnNMSE(2:end);
    ExtraFLOPsFromPreviousPoint(2:end) = diff(pnnnFLOPs);
    nonzeroExtra = ExtraFLOPsFromPreviousPoint(2:end) ~= 0;
    marginal = nan(height(pnnn)-1, 1);
    marginal(nonzeroExtra) = 100*GainFromPreviousPointDb(1 + ...
        find(nonzeroExtra))./ExtraFLOPsFromPreviousPoint(1 + ...
        find(nonzeroExtra));
    MarginalGainDbPer100AdditionalFLOPs(2:end) = marginal;
end
NearOptimal = NMSELossFromBestDb <= ...
    double(selectionConfig.nmseToleranceDb) + 10*eps(abs(bestNMSE));
BeatsGMP = pnnnNMSE <= gmpNMSE + 10*eps(abs(gmpNMSE));
Admissible = NearOptimal & BeatsGMP;
selectedRow = chooseMinimumComplexity( ...
    Admissible, pnnnFLOPs, pnnnParameters);
Selected = false(height(pnnn), 1);
Selected(selectedRow) = true;

diagnostics = pnnn;
diagnostics.NMSELossFromBestDb = NMSELossFromBestDb;
diagnostics.PNNNGainOverGMPDb = PNNNGainOverGMPDb;
diagnostics.FLOPsSavedVsBest = FLOPsSavedVsBest;
diagnostics.FLOPsSavedPercentVsBest = FLOPsSavedPercentVsBest;
diagnostics.ParametersSavedVsBest = ParametersSavedVsBest;
diagnostics.ParametersSavedPercentVsBest = ParametersSavedPercentVsBest;
diagnostics.GainFromPreviousPointDb = GainFromPreviousPointDb;
diagnostics.ExtraFLOPsFromPreviousPoint = ExtraFLOPsFromPreviousPoint;
diagnostics.MarginalGainDbPer100AdditionalFLOPs = ...
    MarginalGainDbPer100AdditionalFLOPs;
diagnostics.NearOptimal = NearOptimal;
diagnostics.BeatsGMP = BeatsGMP;
diagnostics.Admissible = Admissible;
diagnostics.Selected = Selected;

sensitivity = buildSensitivityTable(selectionConfig, pnnnParameters, ...
    pnnnNMSE, pnnnFLOPs, gmpNMSE, bestNMSE, bestRow);
summarySentence = compose([ ...
    'Using a %.2f dB near-optimality tolerance, the lowest-complexity ' ...
    'PNNN configuration within %.2f dB of the best observed PNNN NMSE ' ...
    'and outperforming GMP at the same parameter budget was selected ' ...
    'by the %s. It uses %d active real parameters and %.0f ' ...
    'FLOPs/sample, achieves %.4f dB NMSE, incurs a %.4f dB loss from ' ...
    'the best PNNN point, reduces FLOPs by %.2f%% relative to that best ' ...
    'point, and gains %.4f dB over same-budget Complex GMP-DOMP.'], ...
    double(selectionConfig.nmseToleranceDb), ...
    double(selectionConfig.nmseToleranceDb), criterionName, ...
    pnnnParameters(selectedRow), pnnnFLOPs(selectedRow), ...
    pnnnNMSE(selectedRow), ...
    NMSELossFromBestDb(selectedRow), ...
    FLOPsSavedPercentVsBest(selectedRow), ...
    PNNNGainOverGMPDb(selectedRow));

selection = struct( ...
    'selectedParameters', pnnnParameters(selectedRow), ...
    'criterionName', criterionName, ...
    'nmseToleranceDb', double(selectionConfig.nmseToleranceDb), ...
    'bestPNNNParameters', pnnnParameters(bestRow), ...
    'bestPNNNNMSEdB', bestNMSE, ...
    'selectedPNNNNMSEdB', pnnnNMSE(selectedRow), ...
    'selectedPNNNFLOPs', pnnnFLOPs(selectedRow), ...
    'nmseLossFromBestDb', NMSELossFromBestDb(selectedRow), ...
    'flopsSavedPercentVsBest', ...
        FLOPsSavedPercentVsBest(selectedRow), ...
    'pnnnGainOverGMPDb', PNNNGainOverGMPDb(selectedRow), ...
    'diagnosticsTable', diagnostics, ...
    'sensitivityTable', sensitivity, ...
    'summarySentence', string(summarySentence));
end

function selectedRow = chooseMinimumComplexity(admissible, flops, parameters)
candidateRows = find(admissible);
if isempty(candidateRows)
    error('selectOperatingPoint:NoAdmissiblePoint', ...
        ['No PNNN point is both near-optimal and no worse than its ' ...
        'same-budget Complex GMP pair.']);
end
[~, order] = sortrows([flops(candidateRows), parameters(candidateRows)], ...
    [1 2]);
selectedRow = candidateRows(order(1));
end

function sensitivity = buildSensitivityTable(config, parameters, ...
    nmse, flops, gmpNMSE, bestNMSE, bestRow)
tolerances = double(config.sensitivityTolerancesDb(:));
count = numel(tolerances);
NMSEToleranceDb = tolerances;
SelectedParameters = zeros(count, 1);
SelectedPNNNNMSEdB = zeros(count, 1);
SelectedPNNNFLOPs = zeros(count, 1);
NMSELossFromBestDb = zeros(count, 1);
FLOPsSavedPercentVsBest = zeros(count, 1);
PNNNGainOverGMPDb = zeros(count, 1);
for index = 1:count
    admissible = nmse - bestNMSE <= tolerances(index) + ...
        10*eps(abs(bestNMSE)) & nmse <= gmpNMSE + 10*eps(abs(gmpNMSE));
    row = chooseMinimumComplexity(admissible, flops, parameters);
    SelectedParameters(index) = parameters(row);
    SelectedPNNNNMSEdB(index) = nmse(row);
    SelectedPNNNFLOPs(index) = flops(row);
    NMSELossFromBestDb(index) = nmse(row) - bestNMSE;
    FLOPsSavedPercentVsBest(index) = ...
        100*(flops(bestRow) - flops(row))/flops(bestRow);
    PNNNGainOverGMPDb(index) = gmpNMSE(row) - nmse(row);
end
sensitivity = table(NMSEToleranceDb, SelectedParameters, ...
    SelectedPNNNNMSEdB, SelectedPNNNFLOPs, NMSELossFromBestDb, ...
    FLOPsSavedPercentVsBest, PNNNGainOverGMPDb);
end
