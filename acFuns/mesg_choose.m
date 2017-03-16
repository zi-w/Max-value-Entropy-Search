% Author: Zi Wang
function optimum = mesg_choose(nM, nK, xx, yy, KernelMatrixInv, ...
    guesses, sigma0, sigma, l, xmin, xmax, epsilon)
% This function returns the next evaluation point using MES-G.
% nM is the number of sampled GP hyper-parameter settings.
% nK is the number of sampled maximum values.
% xx, yy are the current observations.
% KernelMatrixInv is the gram maxtrix inverse under different GP
% hyper-parameters.
% guesses are the inferred points to recommend to evaluate.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.
% epsilon is an offset on the sampled max-value.
if nargin <= 11; epsilon = 0.1; end
% Sample a set of random points in the search space.
gridSize = 10000;
d = size(xmin, 1);
Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, 1) ...
    .* rand(gridSize, d);
Xgrid = [ Xgrid ; guesses; xx ];
sx = size(Xgrid, 1);
maxes = zeros(nM, nK);
yvals = 0;
for i = 1:nM
    % Compute the posterior mean/variance predictions for Xgrid.
    [meanVector, varVector] = mean_var(Xgrid, xx, yy, ...
        KernelMatrixInv{i}, l(i,:), sigma(i), sigma0(i));
    % Avoid numerical errors by enforcing variance to be positive.
    varVector(varVector<=sigma0(i)+eps) = sigma0(i)+eps;
    % Obtain the posterior standard deviation.
    sigVector = sqrt(varVector);
    % Define the CDF of the function upper bound.
    probf = @(m0) prod(normcdf((m0 - meanVector)./sigVector));
    % Randomly sample the function upper bounds from a Gumbel distribution
    % that approximates the CDF.
    
    % Find the sample range [left, right].
    
    % Use the fact that the function upper bound is greater than the max of
    % the observations.
    left = max(yy);
    if probf(left) < 0.25
        right = max(meanVector+5*sigVector);
        while (probf(right) < 0.75)
            right = right + right - left;
        end
        mgrid = linspace(left, right, 100);
        
        prob = prod(normcdf((repmat(mgrid,[sx,1]) - repmat(meanVector ...
            ,[1,100]))./repmat(sigVector,[1, 100])),1);
        % Find the median and quartiles.
        med = find_between(0.5, probf, prob, mgrid, 0.01);
        q1 = find_between(0.25, probf, prob, mgrid, 0.01);
        q2 = find_between(0.75, probf, prob, mgrid, 0.01);
        % Approximate the Gumbel parameters alpha and beta.
        beta=(q1-q2)/(log(log(4/3)) - log(log(4)));
        alpha = med+beta*log(log(2));
        assert(beta > 0);
        % Sample from the Gumbel distribution.
        maxes(i,:) = - log( -log(rand(1, nK)) ) .* beta + alpha;
        maxes(i, maxes(i,:) < left + epsilon) = epsilon;
    else
        % In rare cases, the GP shows that with probability at least 0.25, 
        % the function upper bound is smaller than the max of
        % the observations. We manually set the samples maxes to be 
        maxes(i,:) = left + epsilon;
    end
    % Compute the acquisition function values on Xgrid.
    gamma = (repmat(maxes(i,:),[sx 1]) - repmat(meanVector, [1, nK])) ...
        ./ repmat(sigVector, [1, nK]);
    pdfgamma = normpdf(gamma);
    cdfgamma = normcdf(gamma);
    yvals = yvals + sum(gamma.*pdfgamma./(2*cdfgamma) - log(cdfgamma),2);
end

[maxVal, maxIdx] = max(yvals);
start = Xgrid(maxIdx,:);
% Optimize the acquisition function.
acfun = @(x) evaluateMES(x, maxes, xx, yy, KernelMatrixInv, l, sigma, sigma0);
neg_acfun = @(x) negative_wrapper(acfun, x);
[ optimum, fval] = fmincon(neg_acfun, start, [], [], [], [], xmin, xmax, ... 
    [], optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', ...
    'GradObj', 'on'));
fval = -fval;
if fval < maxVal
    disp('Optimization for MES-G failed.')
    optimum = start;
end
