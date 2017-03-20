% Copyright (c) 2017 Zi Wang
function optimum = est_choose(nM, xx, yy, KernelMatrixInv, guesses, ...
    sigma0, sigma, l, xmin, xmax)
% This function returns the next evaluation point using EST.
% nM is the number of sampled GP hyper-parameter settings.
% xx, yy are the current observations.
% KernelMatrixInv is the gram maxtrix inverse under different GP
% hyper-parameters.
% guesses are the inferred points to recommend to evaluate.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.

% Sample a set of random points in the search space.
gridSize = 10000;
d = size(xmin, 1);
Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, ...
    1) .* rand(gridSize, d);
Xgrid = [ Xgrid ; guesses; xx];
maxes = zeros(nM, 1);
yvals = 0;
for i = 1:nM
    % Compute the posterior mean/variance predictions for Xgrid.
    [meanVector, varVector] = mean_var(Xgrid, xx, yy, KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
    % Avoid numerical errors by enforcing variance to be positive.
    varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
    % Obtain the posterior standard deviation.
    sigVector = sqrt(varVector);
    % Compute the expected maximum of the function.
    m = max(yy);
    m0 = m;
    intcnt = 1;
    logprodphi = zeros(100,1);
    while intcnt == 1 || logprodphi(intcnt-1) < 0
        logprodphi(intcnt) = sum(log(normcdf((m0-meanVector)./sigVector, 0, 1)))';
        m = m + (1-exp(logprodphi(intcnt)))*0.05;
        m0 = m0+0.05;
        intcnt = intcnt + 1;
    end
    maxes(i) = m;
    % Compute the EST acquisition function values.
    yvals = yvals + (m - meanVector) ./ sigVector;
end
yvals = yvals / nM;

[minVal, minIdx] = min(yvals);
start = Xgrid(minIdx,:);

% Optimize the acquisition function.
acfun = @(x) evaluateEST(x, xx, yy, KernelMatrixInv, l, sigma, sigma0, maxes);
[ optimum, fval] = fmincon(acfun, start, [], [], [], [], xmin, xmax, [], ...
    optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', 'GradObj', 'on'));
if minVal < fval
    disp('fmincon for EST failed to return a value better than the initialization.')
    [optimum, fval] = fminsearch(acfun, start);
    if minVal < fval
        disp('fminsearch for EST failed to return a value better than the initialization.')
        optimum = start;
    else
        % check if optimum is in range
        flag = 0;
        for i = 1:length(xmin)
            if optimum(i) < xmin(i) || optimum(i) > xmax(i)
                flag = 1;
                break;
            end
        end
        if flag
            disp('fminsearch for EST failed to search within [xmin, xmax].')
            optimum = start;
        end
    end
end