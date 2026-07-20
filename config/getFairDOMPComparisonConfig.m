function cfg = getFairDOMPComparisonConfig(projectRoot)
% getFairDOMPComparisonConfig - Configure the full-signal comparison.
% GMP and PNNN share one deterministic 10% identification domain.

if nargin < 1 || isempty(projectRoot)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
end
projectRoot = char(string(projectRoot));

cfg = struct();
cfg.measurementFile = fullfile(projectRoot, 'measurements', ...
    'experiment20260429T134032_forward_xy.mat');
cfg.mappingMode = 'xy_forward';

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
cfg.gmp.dompOptions = struct('columnTolerance', 1e-12);
cfg.lambdaGrid = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2];
cfg.fixedRidgeLambdas = [1e-3 1e-4 1e-5];
cfg.reducedRealParameterTarget = 200;

cfg.sweep = struct();
cfg.sweep.schemaVersion = 3;
cfg.sweep.coefficientRangeDefinition = "unit_peak_io_unit_column_norm_v1";
cfg.sweep.parameterGrid = 20:10:500;
cfg.sweep.resume = true;
cfg.sweep.candidateBlockSize = 2048;
cfg.sweep.resultsRoot = fullfile(projectRoot, 'results', 'parameter_sweep');

cfg.pnnn = struct();
cfg.pnnn.M = 13;
cfg.pnnn.orders = [1 3 5 7];
cfg.pnnn.featMode = 'full';
cfg.pnnn.removeDC = true;
cfg.pnnn.sparseBaseHiddenNeurons = 12;
cfg.pnnn.nnSeed = 42;

cfg.training = struct();
cfg.training.optimizer = "adam";
cfg.training.miniBatchSize = 1024;
cfg.training.initialLearnRate = 2e-4;
cfg.training.learnRateSchedule = "piecewise";
cfg.training.learnRateDropFactor = 0.95;
cfg.training.shuffle = "every-epoch";
cfg.training.executionEnvironment = "cpu";
cfg.training.verbose = false;
cfg.training.inputDataFormats = "BC";
cfg.training.targetDataFormats = "BC";
cfg.training.historicalTrainFraction = 0.70;
cfg.training.historicalMaxEpochs = 150;
cfg.training.historicalLearnRateDropPeriod = 5;
cfg.training.historicalValidationPatience = 50;

cfg.pruning = struct();
cfg.pruning.historicalFineTuneEpochs = 20;
cfg.pruning.fineTuneInitialLearnRate = cfg.training.initialLearnRate;
cfg.pruning.fineTuneSeedOffset = 100000;

cfg.selection = struct();
cfg.selection.stabilizationWindowParameters = 100;
cfg.selection.stabilizationToleranceDb = 0.20;
cfg.selection.sensitivityWindowsParameters = [80 100 120];
cfg.selection.sensitivityTolerancesDb = [0.15 0.20 0.25];

cfg.paper = struct();
cfg.paper.matlab2tikzSource = fullfile(projectRoot, 'third_party', ...
    'matlab2tikz', 'src');
cfg.paper.latexmkCommand = 'latexmk';

end
