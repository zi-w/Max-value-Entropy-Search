% Author: Zi Wang
function K = compute_gram(xx, hyp, hyp_idx, z)
% This function computes the gram matrix for an add-GP.
all_cat = unique(z);
K = 0;
for i = 1:length(all_cat)
    K = K + computeKmm(xx(:,z==all_cat(i)), hyp.l(hyp_idx, z==all_cat(i))', ...
        hyp.sigma(hyp_idx, all_cat(i)), hyp.sigma0(hyp_idx, all_cat(i)));
end