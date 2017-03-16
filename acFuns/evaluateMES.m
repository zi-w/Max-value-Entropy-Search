% Author: Zi Wang
function [f, g] = evaluateMES(x, maxes, xx, yy, KernelMatrixInv, ...
    l, sigma, sigma0)
% This function computes the acquisition function value f and gradient g at 
% the queried point x using MES given sampled function values maxes, and
% observations xx (size T x d), yy (size T x 1).
% See also: mean_var.m
if nargout == 2
    % Initialize f, g
    f = 0;
    g = 0;
    sx = size(x,1);
    % Notice that if gradient is queried, sx == 1.
    dx = size(x,2);
    nK = size(maxes,2);
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions and gradients.
        [meanVector, varVector, meangrad, vargrad] = mean_var(x, xx, ...
            yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute gamma.
        meanVector = repmat(meanVector, [1, nK]);
        sigVector = repmat(sigVector, [1, nK]);
        meangrad = repmat(meangrad, [1, nK]);
        vargrad = repmat(vargrad, [1, nK]);
        gamma = (repmat(maxes(i,:),[sx 1]) - meanVector) ./ sigVector;
        % Compute the acquisition function of MES.
        pdfgamma = realnormpdf(gamma);
        cdfgamma = normcdf(gamma);
        f = f + sum(gamma.*pdfgamma./(2*cdfgamma) - log(cdfgamma),2);
        % Compute the gradient of the acquisition function of MES.
        tmp = pdfgamma./cdfgamma;
        gamma = repmat(gamma, [dx, 1]);
        tmp = repmat(tmp, [dx, 1]);
        sigVector = repmat(sigVector, [dx, 1]);
        g = g + sum(((-1-gamma.^2).*tmp/2 - gamma.*tmp.^2/2) .* ...
            (-(0.5*gamma.*vargrad./sigVector + meangrad)./sigVector ), 2);
        % Another option to compute g is
        % g = g + ((-1-gamma.^2)*cdfgamma - gamma*pdfgamma)*pdfgamma * ...
        %  (-(0.5*gamma*vargrad./sigVector + meangrad)./sigVector )/(2*cdfgamma^2);  
    end
    % Average f and g
    f = f / size(KernelMatrixInv, 2) / size(maxes,2);
    g = g / size(KernelMatrixInv, 2) / size(maxes,2);

else
    f = 0;
    sx = size(x,1);
    nK = size(maxes,2);
    for i = 1 :  size(KernelMatrixInv, 2)
        % Compute the posterior mean/variance predictions.
        [meanVector, varVector] = mean_var(x, xx, yy, ...
            KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
        % Avoid numerical errors by enforcing variance to be positive.
        varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
        % Obtain the posterior standard deviation.
        sigVector = sqrt(varVector);
        % Compute gamma.
        meanVector = repmat(meanVector, [1, nK]);
        sigVector = repmat(sigVector, [1, nK]);
        gamma = (repmat(maxes(i,:),[sx 1]) - meanVector) ./ sigVector;
        % Compute the acquisition function of MES.
        pdfgamma = realnormpdf(gamma);
        cdfgamma = normcdf(gamma);
        f = f + sum(gamma.*pdfgamma./(2*cdfgamma) - log(cdfgamma),2);
    end
    % Average f.
    f = f / size(KernelMatrixInv, 2) / size(maxes,2);
end