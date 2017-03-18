% Author: Zi Wang
function nlz = compute_nlz_wrap(xx, yy, params, n_partition, z)
% This function returns the negative log data liklihood for add-GP.

% params(1:n_partition) defines l in log scale
% assume l is the same for all dimensions in each partition
params = params';
l = exp(params(1, z));
% params(n_partition+1:2*n_partition) defines sigma in log scale
sigma = exp(params(1, n_partition+1:2*n_partition));
% params(2*n_partition+1:3*n_partition) defines sigma0 in log scale
sigma0 = exp(params(1, 2*n_partition+1:3*n_partition));
hyp.l = l;
hyp.sigma = sigma; 
hyp.sigma0 = sigma0;

nlz = compute_nlz(xx, yy, hyp, z);