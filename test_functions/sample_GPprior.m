% Copyright (c) 2017 Zi Wang
function f = sample_GPprior(dx)
nFeatures = 10000;
l = ones(1, dx)*256.;
sigma = 5;
sigma0 = 0.0001;
W = randn(nFeatures, dx) .* repmat(sqrt(l), nFeatures, 1);
b = 2 * pi * rand(nFeatures, 1);
theta = randn(nFeatures, 1);
f = @(x) (theta' * sqrt(2 * sigma / nFeatures) * cos(W * x' + repmat(b, 1, size(x, 1))))';
