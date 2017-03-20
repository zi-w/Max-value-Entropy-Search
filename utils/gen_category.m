% Author: Zi Wang
function z = gen_category(dx)
% Generate synthetic category,
% with at most 3 dimensional substructure
z = zeros(1,dx);
cnt = 0;
groupidx = 1;
while cnt < dx
    randgroup = randi([1,3]);
    if rand < 0.5 && randgroup == 1
        randgroup = 2;
    end
    z(cnt+1:cnt+randgroup) = groupidx;
    groupidx = groupidx + 1;
    cnt = cnt + randgroup;
end
z = z(1:dx);