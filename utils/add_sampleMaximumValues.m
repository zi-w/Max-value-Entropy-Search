% Author: Zi Wang
% See also: sampleMaximumValues.m
function [ samples] = add_sampleMaximumValues(nK, xx, yy, sigma0, ...
    sigma, l, xmin, xmax, nFeatures, cz, all_cat)
% This function returns sampled maximum values for each component of the
% posterior add-GP conditioned on current obervations. We construct random
% features and optimize function components drawn from the posterior GP.
% nK is the number of sampled maximum values.
% xx, yy are the current observations.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.
% nFeatures is the number of random features sampled to approximate the
% GP.
% cz is the category assignment vector of all input dimensions.
% all_cat is the set of unique categories in cz.
nZ = length(all_cat);
nX = size(xx, 1);
assert(nX > 0 && nX < nZ*nFeatures)
samples = zeros(nZ, nK)*-1e10;

% Draw weights for the random features of each function component.
W = cell(1,nZ);
b = cell(1,nZ);
% Random features for xx.
Z = zeros(nFeatures * nZ, nX);
for k = 1:nZ
    coords = (cz==all_cat(k));
    d = sum( coords );
    W{k} = randn(nFeatures, d) .* repmat(sqrt(l(1,coords)), nFeatures, 1);
    b{k} = 2 * pi * rand(nFeatures, 1);
    % Compute the features for xx of the k-th function component.
    Z((k-1)*nFeatures+1:k*nFeatures,:) = sqrt(2 * sigma(1,all_cat(k)) ...
        / (nZ*nFeatures)) * cos(W{k} * xx(:,coords)' + repmat(b{k}, 1, nX));
end

% Compute the mean and covariance matrix of the coefficient theta.
sigma0_all = sum(sigma0(1,all_cat));
if (nX < nZ*nFeatures)
    % We adopt the formula $theta \sim \N(Z(Z'Z + \sigma^2 I)^{-1} y,
    % I-Z(Z'Z + \sigma^2 I)Z')$.
    Sigma = Z' * Z + sigma0_all * eye(size(xx, 1));
    mu = Z*chol2invchol(Sigma)*yy;
    [U, D] = eig(Sigma);
    D = diag(D);
    R = (sqrt(D) .* (sqrt(D) + sqrt(sigma0_all))).^-1;
else
    % $theta \sim \N((ZZ'/\sigma^2 + I)^{-1} Z y / \sigma^2,
    % (ZZ'/\sigma^2 + I)^{-1})$.
    Sigma = chol2invchol(Z*Z' / sigma0_all + eye(nZ*nFeatures));
    mu = Sigma * Z * yy / sigma0_all;
end
for j = 1:nK
    % Sample the coefficient theta.
    noise = randn(nZ*nFeatures, 1);
    if (nX < nZ*nFeatures)
        theta = noise - (Z * (U * (R .* (U' * (Z' * noise))))) + mu;
    else
        theta = mu + noise * chol(Sigma);
    end
    for k = 1:nZ
        coords = (cz==all_cat(k));
        d = sum( coords );
        theta_sub = theta((k-1)*nFeatures+1:k*nFeatures);
        xmin_sub = xmin(coords);
        xmax_sub = xmax(coords);
        % Obtain a function component sampled from the posterior GP.
        targetVector = @(x) (theta_sub' * sqrt(2 * sigma(1,all_cat(k)) / (nZ*nFeatures)) * cos(W{k} * x' + repmat(b{k}, 1, size(x, 1))))';
        targetVectorGradient = @(x) theta_sub' * -sqrt(2 * sigma(1,all_cat(k)) / (nZ*nFeatures)) * (repmat(sin(W{k} * x' + b{k}), 1, d) .* W{k});
        target = @(x) wrap_target(targetVector, targetVectorGradient, x);
        % Optimize the function component.
        [~, sample]= globalMaximization(target, xmin_sub, xmax_sub, xx(:,coords));
        samples(k, j) = sample;
    end
end

end


function [f,g] = wrap_target(tf, tg, x)
f = tf(x);
if nargout > 1
    g = tg(x);
end
end