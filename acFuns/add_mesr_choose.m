% Copyright (c) 2017 Zi Wang
function optimum = add_mesr_choose(maxes, xx, yy, KernelMatrixInv, sigma0, sigma, l, xmin, xmax)
% This function returns the next evaluation point using add-MES-R.
% maxes is a nZ x nK matrix with sampled max-values of size nK for each 
% function component of the add-GP with nZ components. 
% xx, yy are the current observations.
% KernelMatrixInv is the gram maxtrix inverse under different GP
% hyper-parameters.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.

% Define the acquisition function (and gradient) of MES.
acfun = @(x) evaluateMES(x, maxes, xx, yy, KernelMatrixInv, l, sigma, sigma0);
% We optimize globally the cost function
% Optimize the acquisition function.
optimum = globalMaximization(acfun, xmin, xmax, xx);