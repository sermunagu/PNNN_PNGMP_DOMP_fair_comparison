function [X_in, r_vec] = buildPhaseNormInput(x, M, orders, featMode)
% buildPhaseNormInput - Build phase-normalized PNNN input features.
%
% This function rotates each sample to the local input phase, builds the
% configured memory/nonlinear feature vector, and returns the rotations used
% later to reconstruct the complex model output.
%
% Inputs:
%   x - Modeled-block input signal under the local X/Y convention.
%   M, orders, featMode - Memory depth, nonlinear orders, and feature mode.
%
% Outputs:
%   X_in - D x N phase-normalized input feature matrix.
%   r_vec - 1 x N phase-rotation vector.

    if nargin < 4 || isempty(featMode)
        featMode = "full";
    end

    x = x(:).';

    validateInputSignal(x, M);
    orders = validateOrders(orders);
    featMode = string(featMode);

    N = numel(x);

    switch featMode
        case "full"
            D = 2*(M+1) + (M+1)*numel(orders);
        case "pruned"
            D = 1 + 2*M + M*numel(orders);
        otherwise
            error("featMode debe ser 'full' o 'pruned'.");
    end

    X_in = zeros(D, N);
    r_vec = complex(zeros(1, N));

    for k = 1:N
        if abs(x(k)) == 0
            r = 1 + 0j;
        else
            r = conj(x(k)) / abs(x(k));
        end

        taps = periodicTapIndices(k, M, N);
        x_k = x(taps);
        X_k = r * x_k;

        switch featMode
            case "full"
                re_X_k = real(X_k);
                im_X_k = imag(X_k);
                A_k = abs(x_k);

                env = zeros(1, (M+1)*numel(orders));
                for ip = 1:numel(orders)
                    cols = (ip-1)*(M+1) + (1:M+1);
                    env(cols) = A_k.^orders(ip);
                end

                phi_k = [re_X_k, im_X_k, env];

            case "pruned"
                re_main = real(X_k(1));
                re_mem = real(X_k(2:end));
                im_mem = imag(X_k(2:end));
                A_mem = abs(x_k(2:end));

                env = zeros(1, M*numel(orders));
                for ip = 1:numel(orders)
                    cols = (ip-1)*M + (1:M);
                    env(cols) = A_mem.^orders(ip);
                end

                phi_k = [re_main, re_mem, im_mem, env];
        end

        X_in(:, k) = phi_k(:);
        r_vec(k) = r;
    end
end

function taps = periodicTapIndices(k, M, N)
    taps = mod((k:-1:k-M) - 1, N) + 1;
end

function validateInputSignal(x, M)
    if isempty(x)
        error('x no puede estar vacío.');
    end
    if any(~isfinite(x))
        error('x contiene NaN/Inf.');
    end
    if ~isscalar(M) || M < 0 || M ~= floor(M)
        error('M debe ser un entero no negativo.');
    end
    if numel(x) <= M
        error('La señal debe tener más muestras que M. N=%d, M=%d.', numel(x), M);
    end
end

function orders = validateOrders(orders)
    orders = orders(:).';
    if isempty(orders)
        error('orders no puede estar vacío.');
    end
    if any(~isfinite(orders)) || any(orders < 1) || any(orders ~= floor(orders))
        error('orders debe contener enteros positivos finitos.');
    end
end
