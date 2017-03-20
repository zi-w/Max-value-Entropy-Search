% Copyright (c) 2017 Zi Wang
function K = computeAddKnm(xx, xxp, hyp, hyp_idx, z)
% Compute the additive kernel Knm 
all_cat = unique(z);
K = 0;
for i = 1:length(all_cat)
    K = K + computeKnm(xx(:,z==all_cat(i)), xxp(:,z==all_cat(i)), ...
        hyp.l(hyp_idx, z==all_cat(i))', hyp.sigma(hyp_idx, all_cat(i)));
end