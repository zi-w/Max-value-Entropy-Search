function [optimum, fval] = pi_choose(xx, yy, KernelMatrixInv, sigma0, sigma, l, xmin, xmax, m0)
% We sample from the global minimum

    % We define the cost function to be optimized
    %m0 = max(yy) + 0.1;
    target_gradient = @(x) evaluatePIAc(x, xx, yy, KernelMatrixInv, l, sigma, sigma0, m0);
    target = @(x) evaluatePIAc_target(x, xx, yy, KernelMatrixInv, l, sigma, sigma0, m0);
    % We optimize globally the cost function

    [optimum, fval] = globalMaximization2(target, target_gradient, xmin, xmax);