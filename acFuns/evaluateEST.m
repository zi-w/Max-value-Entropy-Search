function [f, g] = evaluateEST(x, xx, yy, KernelMatrixInv, l, sigma, ...
    sigma0, maxes)
% This function computes the acquisition function value f and gradient g at
% the queried point x using EST given sampled function values maxes, and
% observations xx (size T x d), yy (size T x 1).

f = 0;
g = 0;
for i = 1 :  size(KernelMatrixInv, 2)
    % Compute the posterior mean/variance predictions and gradients.
    [meanVector, varVector, meangrad, vargrad] = mean_var(x, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
    % Avoid numerical errors by enforcing variance to be positive.
    varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
    % Obtain the posterior standard deviation.
    sigVector = sqrt(varVector);
    % Compute the EST acquisition function values f and gradients g.
    gamma = (maxes(i) - meanVector)./sigVector;
    f = f + gamma;
    g = g + (-(0.5*gamma*vargrad./sigVector + meangrad)./sigVector );
end
f = f / size(KernelMatrixInv, 2);
g = g / size(KernelMatrixInv, 2);
