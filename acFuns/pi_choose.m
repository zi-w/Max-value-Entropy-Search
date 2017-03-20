% Copyright (c) 2017 Zi Wang
function [optimum, fval] = pi_choose(xx, yy, KernelMatrixInv, guesses, ...
    sigma0, sigma, l, xmin, xmax)
% This function returns the next evaluation point using PI.
% xx, yy are the current observations.
% KernelMatrixInv is the gram maxtrix inverse under different GP
% hyper-parameters.
% guesses are the inferred points to recommend to evaluate.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.

% Define the acquisition function (and gradient) of EI.
acfun = @(x) evaluatePI(x, xx, yy, KernelMatrixInv, l, sigma, sigma0);
% Optimize the acquisition function.
[optimum, fval] = globalMaximization(acfun, xmin, xmax, [guesses;xx]);