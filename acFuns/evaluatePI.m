% Copyright (c) 2017 Zi Wang
function [f, g] = evaluatePI(x, xx, yy, KernelMatrixInv, l, sigma, sigma0)
% This function computes the acquisition function value f and gradient g at
% the queried point x using PI given sampled function values maxes, and
% observations Xsamples (size T x d), Ysamples (size T x 1).

if nargout == 2
    f = 0;
    g = 0;
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions and gradients.
        [meanVector, varVector, meangrad, vargrad] = mean_var(x, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the PI acquisition function values f and gradients g.
        gamma = (max(yy) + sqrt(sigma0(i)) - meanVector)./sigVector;
        f = f + gamma;
        g = g + (-(0.5*gamma*vargrad./sigVector + meangrad)./sigVector );
    end
    % Set f, g to negative to maximize later.
    f = - f / size(KernelMatrixInv, 2);
    g = - g / size(KernelMatrixInv, 2);
else
    f = 0;
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions and gradients.
        [meanVector, varVector] = mean_var(x, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the PI acquisition function values f.
        f = f + (max(yy) + sqrt(sigma0(i)) - meanVector)./sigVector;
    end
    % Set f to negative to maximize later.
    f = - f / size(KernelMatrixInv, 2);
end