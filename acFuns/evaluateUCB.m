% Copyright (c) 2017 Zi Wang
function [f, g] = evaluateUCB(x, xx, yy, KernelMatrixInv, l, sigma, ...
    sigma0, alpha, beta)
% This function computes the acquisition function value f and gradient g at
% the queried point x using GP-UCB given sampled function values maxes, and
% observations xx (size T x d), yy (size T x 1).
if nargout == 2
    f = 0;
    g = 0;
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions and gradients.
        [meanVector, varVector, meangrad, vargrad] = mean_var(x, xx, yy, ...
            KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the UCB acquisition function values f and gradients g.
        f = f + alpha * meanVector + beta * sigVector;
        g = g + alpha * meangrad + beta * 0.5 * vargrad./sigVector;
    end
    f = f / size(KernelMatrixInv, 2);
    g = g / size(KernelMatrixInv, 2);
else
    f = 0;
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions.
        [meanVector, varVector] = mean_var(x, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=0) = eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the UCB acquisition function values f.
        f = f + alpha * meanVector + beta * sigVector;
    end
    f = f / size(KernelMatrixInv, 2);
end