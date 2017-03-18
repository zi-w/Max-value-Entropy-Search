% Author: Zi Wang
function z = sample_categorical(p, dx)
% This function samples from a categorical distribution 
rd = rand(1,dx);
cur_p = 0;
z = zeros(1,dx);
for i = 1:length(p)
    cur_p = cur_p + p(i);
    z = z + (z==0 & rd <= cur_p)*i;
end
