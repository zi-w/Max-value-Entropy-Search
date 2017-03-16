function [f, g] = evaluateEI(x, xx, yy, KernelMatrixInv, l, sigma, sigma0)
% This function computes the acquisition function value f and gradient g at
% the queried point x using EI given sampled function values maxes, and
% observations Xsamples (size T x d), Ysamples (size T x 1).
if nargout == 2
    f = 0;
    g = 0;
    m0 = max(yy);
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions and gradients.
        [meanVector, varVector, meangrad, vargrad] = mean_var(x, xx, yy, ...
            KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the EI acquisition function values f and gradients g.
        gamma = (m0 - meanVector)./sigVector;
        pdfgamma =normpdf(gamma);
        cdfgamma = normcdf(gamma);
        temp = (pdfgamma - gamma .* (1-cdfgamma));
        f = f + sigVector.*temp;
        g = g + 0.5*temp.*vargrad ./ sigVector + sigVector .* (cdfgamma-1) ...
            .*(-(0.5*gamma*vargrad./sigVector + meangrad)./sigVector );
    end
    f = f / size(KernelMatrixInv, 2);
    g = g / size(KernelMatrixInv, 2);
else
    f = 0;
    m0 = max(Ysamples);
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions.
        [meanVector, varVector] = mean_var(x, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute the EI acquisition function values f.
        gamma = (m0 - meanVector)./sigVector;
        pdfgamma = normpdf(gamma);
        cdfgamma = normcdf(gamma);
        temp = (pdfgamma - gamma .* (1-cdfgamma));
        f = f + sigVector.*temp;
    end
    f = f / size(KernelMatrixInv, 2);
end