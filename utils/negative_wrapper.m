% Author: Zi Wang
function [neg_f, neg_g] = negative_wrapper(f, x)
% Wrap the function f to its negative values
[neg_f, neg_g] = f(x);
neg_f = -neg_f;
neg_g = -neg_g;
