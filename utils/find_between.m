% Author: Zi Wang
function res = find_between(val, func, funcvals, mgrid, thres)
% This function finds res\in\R such that func(res) = val via binary search.
% mgrid is a vector such that min(mgrid) < res < max(mgrid).
% funcvals = func(mgrid).
% thres threshold for the precision of res.
[~,t2]=min(abs(funcvals-val));
if (abs(funcvals(t2) - val) < thres)
    res = mgrid(t2);
    return;
end
assert(funcvals(1) < val && funcvals(end) > val);
if funcvals(t2) > val
    left = mgrid(t2-1);
    right = mgrid(t2);
else
    left = mgrid(t2);
    right = mgrid(t2+1);
end
mid = (left + right)/2;
midval = func(mid);
while abs(midval - val) > thres
    if midval > val
        right = mid;
    else
        left = mid;
    end
    mid = (left + right) / 2;
    midval = func(mid);
end
res = mid;