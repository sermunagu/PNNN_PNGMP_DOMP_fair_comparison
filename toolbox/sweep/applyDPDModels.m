function [dpd, metadata] = applyDPDModels(bundle, u, outputFile)
% Infer a complete periodic record. Never fit, rescale, clip or restore DC.
% Input must have the training sample rate and physical amplitude convention.
if ischar(bundle) || isstring(bundle)
    saved = load(bundle, 'bundle');
    bundle = saved.bundle;
end
if nargin < 3, outputFile = ''; end
if ~isempty(outputFile) && isfile(outputFile)
    error('applyDPDModels:WouldOverwrite', 'File already exists: %s', outputFile);
end
validateattributes(u, {'double','single'}, {'vector','nonempty','finite'});
x = double(u(:));
inputMean = mean(x);
if bundle.config.pnnn.removeDC, x = x - inputMean; end
rows = (1:numel(x)).';
labels = ["Complex_GMP_DOMP", "PN_IQ_GMP", "PNNN"];
dpd = repmat(struct('yvalmod', [], 'modeltype', ''), 1, 3);
for k = 1:3
    model = bundle.models{k};
    if model.kind == "pnnn"
        f = model.featureConfig;
        [features, ~, rotation] = buildPhaseNormDataset( ...
            x, zeros(size(x)), f.M, f.orders, f.featMode);
        stats = model.normalization;
        features = (features.' - stats.muX) ./ stats.sigmaX;
        prediction = predictPhaseNorm(model.network, features, stats, rotation);
    else
        manager = struct('regPopulation', model.regPopulation);
        support = (1:numel(model.regPopulation)).';
        prediction = complex(zeros(size(x)));
        for first = 1:bundle.config.gmp.blockSize:numel(x)
            local = first:min(first + bundle.config.gmp.blockSize - 1, numel(x));
            if model.kind == "complex"
                U = buildGMPRegressorRows(x, rows(local), manager, support);
                prediction(local) = U * model.coefficients;
            else
                rotation = complex(ones(numel(local), 1));
                nonzero = abs(x(local)) ~= 0;
                values = x(local);
                rotation(nonzero) = conj(values(nonzero)) ./ abs(values(nonzero));
                U = buildPNIQFeatures(x, rows(local), rotation, ...
                    manager, support, model.descriptors);
                U = U(:, model.selectedColumns);
                prediction(local) = conj(rotation) .* complex( ...
                    U * model.coefficientsI, U * model.coefficientsQ);
            end
        end
    end
    assert(numel(prediction) == numel(x) && all(isfinite(prediction)), ...
        'Invalid DPD prediction.');
    dpd(k).yvalmod = complex(prediction(:));
    dpd(k).modeltype = char(labels(k) + "_P" + bundle.P);
end
metadata = rmfield(bundle, 'models');
metadata.inferenceInputMean = inputMean;
metadata.inferenceSampleCount = numel(x);
metadata.outputDomain = "DC-removed ILC target; no output DC restoration";
if ~bundle.config.pnnn.removeDC, metadata.outputDomain = "ILC target"; end
metadata.rfValidationPerformed = false;
if ~isempty(outputFile)
    save(outputFile, 'dpd', 'metadata', '-v7.3');
end
end
