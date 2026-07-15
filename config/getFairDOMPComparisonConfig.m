function cfg = getFairDOMPComparisonConfig(projectRoot)
% getFairDOMPComparisonConfig - Configure the full-signal comparison.
% GMP and PNNN share one deterministic 10% identification domain.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end
projectRoot = char(string(projectRoot));

cfg = struct();
cfg.projectRoot = projectRoot;
cfg.measurementName = 'experiment20260429T134032_xy';
cfg.measurementFile = fullfile(projectRoot, 'measurements', ...
    [cfg.measurementName '.mat']);
cfg.mappingMode = 'xy_forward';
cfg.resultsRoot = fullfile(projectRoot, 'results', ...
    'full_signal_domp_comparison');
cfg.historicalDisjointResultDirectory = fullfile(projectRoot, 'results', ...
    'fair_domp_comparison', '20260714_013938');

cfg.identificationFraction = 0.10;
cfg.identificationSeed = 1004;
cfg.internalTrainFraction = 0.85;
cfg.internalSplitSeed = 42;
cfg.amplitudeBinCount = 10;

cfg.gmp = struct();
cfg.gmp.Qpmax = 50;
cfg.gmp.Qnmax = 50;
cfg.gmp.Pmax = 13;
cfg.gmp.maxPopulation = 100;
cfg.gmp.blockSize = 8192;
cfg.gmp.verbosity = 0;
cfg.gmp.selectionMethod = 'DOMP';
cfg.gmp.dompOptions = struct( ...
    'columnTolerance', 1e-12, ...
    'correlationTolerance', 1e-14, ...
    'residualTolerance', 1e-12, ...
    'lsTolerance', []);
cfg.lambdaGrid = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2];
cfg.reducedRealParameterTarget = 200;

cfg.pnnn = struct();
cfg.pnnn.M = 13;
cfg.pnnn.orders = [1 3 5 7];
cfg.pnnn.featMode = 'full';
cfg.pnnn.actType = 'sigmoid';
cfg.pnnn.removeDC = true;
cfg.pnnn.parameterMatchedTargetModel = ...
    'Independent PN-IQ full';
cfg.pnnn.denseControlHiddenNeurons = 4;
cfg.pnnn.sparseBaseHiddenNeurons = 12;
cfg.pnnn.trainHistoricalN25 = false;
cfg.pnnn.nnSeeds = 42;

cfg.training = struct();
cfg.training.optimizer = "adam";
cfg.training.maxEpochs = [];
cfg.training.miniBatchSize = 1024;
cfg.training.initialLearnRate = 2e-4;
cfg.training.learnRateSchedule = "piecewise";
cfg.training.learnRateDropPeriod = [];
cfg.training.learnRateDropFactor = 0.95;
cfg.training.validationPatience = [];
cfg.training.shuffle = "every-epoch";
cfg.training.executionEnvironment = "cpu";
cfg.training.trainingPlots = "none";
cfg.training.verbose = false;
cfg.training.inputDataFormats = "BC";
cfg.training.targetDataFormats = "BC";
cfg.training.historicalTrainFraction = 0.70;
cfg.training.historicalMaxEpochs = 150;
cfg.training.historicalLearnRateDropPeriod = 5;
cfg.training.historicalValidationPatience = 50;

cfg.pruning = struct();
cfg.pruning.enabled = true;
cfg.pruning.sparsity = 0;
cfg.pruning.targetMode = "activeTrainableParams";
cfg.pruning.targetActiveTrainableParams = [];
cfg.pruning.scope = "global";
cfg.pruning.includeBiases = false;
cfg.pruning.structureMode = "unstructured";
cfg.pruning.structuredRanking = "magnitude";
cfg.pruning.structuredTargetPolicy = "closestNotAbove";
cfg.pruning.hybridExactTarget = false;
cfg.pruning.fineTuneEnabled = true;
cfg.pruning.fineTuneEpochs = [];
cfg.pruning.historicalFineTuneEpochs = 20;
cfg.pruning.fineTuneInitialLearnRate = ...
    cfg.training.initialLearnRate;
cfg.pruning.fineTuneLearnRateDropPeriod = [];
cfg.pruning.fineTuneSeedOffset = 100000;
cfg.pruning.freezePruned = true;

cfg.warmStart = struct();
cfg.warmStart.enabled = false;
cfg.warmStart.useLatestDeploy = false;

cfg.equivalence.relativePredictionTolerance = 1e-10;
cfg.equivalence.nmseDifferenceToleranceDb = 1e-9;
cfg.report.compilePDF = false;
end
