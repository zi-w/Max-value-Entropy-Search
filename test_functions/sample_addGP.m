% Copyright (c) 2017 Zi Wang
function f = sample_addGP(dx, n_partition, xmin, xmax)
% Generate synthetic test functions for add gp opt.
z_o = gen_category(dx);
fixhyp.l = ones(1,dx)*50;
fixhyp.sigma = ones(1, n_partition)*5;
fixhyp.sigma0 = 0.0001*ones(1, n_partition);



n_x = 1000;
xfdata = repmat(xmin', n_x, 1) + repmat(xmax' - xmin', n_x, 1) .* rand(n_x, dx);
K = compute_gram(xfdata, fixhyp, 1, z_o);
yfdata = chol(K)'*randn(n_x,1) + randn(n_x,1)*1;%mvnrnd(zeros(n_x,1), K)';

KernelMatrix = K;
KernelMatrixInv = chol2invchol(KernelMatrix);
f = @(x) computeAddKnm(x, xfdata, fixhyp, 1, z_o) * KernelMatrixInv * yfdata;

% Generate initialization
% n_x = 50;
% initx = repmat(xmin', n_x, 1) + repmat(xmax' - xmin', n_x, 1) .* rand(n_x, dx);
% inity = f(initx);

% Save the function
% save([dirname num2str(i) '.mat'], 'xmin', 'xmax', 'initx', 'inity', 'f', 'fixhyp', 'z_o', 'xfdata', 'yfdata')

