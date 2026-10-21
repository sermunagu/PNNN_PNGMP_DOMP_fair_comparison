function diagnostics = diagnosePNIQAblation(ablationFile)
% Read-only selected-matrix audit: actual LS scaling, refit, directed sensitivity.
saved = load(ablationFile,'result'); result = saved.result;
data = load(fullfile(result.sourceDirectory,'experiment.mat'),'x','y','split','cfg');
x = data.x; y = data.y; id = data.split.identificationIndices;
diagnostics = table();
for k = 1:numel(result.models)
    model = result.models{k}; manager = struct('regPopulation',model.regPopulation);
    sources = (1:numel(model.regPopulation)).'; rotation = ones(numel(id),1);
    if model.kind == "pniq"
        nonzero = abs(x(id)) ~= 0;
        rotation(nonzero) = conj(x(id(nonzero)))./abs(x(id(nonzero)));
        A = buildPNIQFeatures(x,id,rotation,manager,sources,model.descriptors);
    else
        A = buildGMPRegressorRows(x,id,manager,sources);
        if model.kind == "iq", A = [real(A),imag(A)]; end
    end
    if model.kind ~= "complex", A = A(:,model.selectedColumns); end
    peaks = max(abs(A),[],1); B = A./peaks;
    target = rotation.*y(id); outputPeak = max(abs(target)); target = target/outputPeak;
    [~,S,V] = svd(B,'econ'); singular = diag(S);
    tolerance = max(size(B))*eps(singular(1));
    % A perturbation of relative size 1e-12 along the weakest left singular direction.
    direction = B*V(:,end); direction = direction/norm(direction);
    perturbed = target + 1e-12*norm(target)*direction;
    if model.kind == "complex"
        c = lsqminnorm(B,target,tolerance);
        cp = lsqminnorm(B,perturbed,tolerance);
        fitted = model; fitted.coefficients = c*outputPeak./peaks.';
        changed = model; changed.coefficients = cp*outputPeak./peaks.';
    else
        c = complex(lsqminnorm(B,real(target),tolerance),lsqminnorm(B,imag(target),tolerance));
        cp = complex(lsqminnorm(B,real(perturbed),tolerance),lsqminnorm(B,imag(perturbed),tolerance));
        fitted = model; fitted.coefficientsI = real(c)*outputPeak./peaks.';
        fitted.coefficientsQ = imag(c)*outputPeak./peaks.';
        changed = model; changed.coefficientsI = real(cp)*outputPeak./peaks.';
        changed.coefficientsQ = imag(cp)*outputPeak./peaks.';
    end
    replay = predictPNIQAblation(model,x,data.cfg.gmp.blockSize);
    assert(isequal(replay,result.predictions(:,k)));
    refitted = predictPNIQAblation(fitted,x,data.cfg.gmp.blockSize);
    prediction = predictPNIQAblation(changed,x,data.cfg.gmp.blockSize);
    evaluation = data.split.evaluationIndices;
    deltaNMSE = nmseComplexDb(y(evaluation),prediction(evaluation)) - ...
        nmseComplexDb(y(evaluation),refitted(evaluation));
    diagnostics = [diagnostics; table(result.P,result.metrics.Model(k),size(B,2), ...
        nnz(singular>tolerance),singular(1),singular(end),singular(1)/singular(end), ...
        singular(end)/singular(1)<1e-8,max(vecnorm(B))/min(vecnorm(B)), ...
        norm(refitted-replay)/norm(replay),norm(cp-c)/norm(c), ...
        norm(prediction-refitted)/norm(refitted),deltaNMSE, ...
        'VariableNames',{'P','Model','Columns','Rank','SigmaMax','SigmaMin','Condition2', ...
        'NearDependence','ColumnNormRatio','RefitRelativePredictionError', ...
        'PerturbedRelativeCoefficientChange','PerturbedRelativePredictionChange', ...
        'PerturbedEvaluationNMSEChangeDB'})]; %#ok<AGROW>
    fprintf('P%d %s: rank %d/%d, kappa %.3g, sensitivity %.3g dB\n', ...
        result.P,result.metrics.Model(k),nnz(singular>tolerance),size(B,2), ...
        singular(1)/singular(end),deltaNMSE);
end
end
