% Author: Zi Wang
function nlz = compute_nlz(xx, yy, hyp, z)
% Computes the negative log marginal likelihood for add-GP
nlz = 0;
for idx = 1:size(hyp.l,1)
K = compute_gram(xx, hyp, idx, z);
L = chol(K);
alpha = solve_chol(L,yy);
nlz = nlz + yy'*alpha/2 + sum(log(diag(L)));
end