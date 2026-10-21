function [prediction,auxiliary] = predictPNIQAblation(model,x,blockSize)
% Experimental replay on an already DC-preprocessed periodic input, without fitting.
% Auxiliary columns report frozen conjugate / envelope-only contributions.
manager = struct('regPopulation',model.regPopulation);
support = (1:numel(model.regPopulation)).';
degree = arrayfun(@(r) numel(r.X)-numel(r.Xconj),model.regPopulation);
if model.kind ~= "complex", degree = degree(mod(model.selectedColumns-1,numel(support))+1); end
prediction = complex(zeros(size(x))); auxiliary = complex(zeros(numel(x),2));
powers = [-1 0];
for first = 1:blockSize:numel(x)
    rows = (first:min(first+blockSize-1,numel(x))).';
    rotation = ones(numel(rows),1);
    if model.kind == "pniq"
        nonzero = abs(x(rows)) ~= 0;
        rotation(nonzero) = conj(x(rows(nonzero))) ./ abs(x(rows(nonzero)));
        A = buildPNIQFeatures(x,rows,rotation,manager,support,model.descriptors);
    else
        A = buildGMPRegressorRows(x,rows,manager,support);
        if model.kind == "iq", A = [real(A),imag(A)]; end
    end
    if model.kind == "complex"
        prediction(rows) = A*model.coefficients;
    else
        A = A(:,model.selectedColumns);
        prediction(rows) = conj(rotation).*complex(A*model.coefficientsI,A*model.coefficientsQ);
    end
    for a = 1:2
        keep = degree == powers(a);
        if model.kind == "complex"
            auxiliary(rows,a) = A(:,keep)*model.coefficients(keep);
        else
            auxiliary(rows,a) = conj(rotation).*complex( ...
                A(:,keep)*model.coefficientsI(keep),A(:,keep)*model.coefficientsQ(keep));
        end
    end
end
end
