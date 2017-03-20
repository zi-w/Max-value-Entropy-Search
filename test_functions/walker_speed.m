% Copyright (c) 2017 Zi Wang
function ret = walker_speed(hyperparams)
% Wrapper of walker_main.
% hyperparams(1) \in [1,10]
% hyperparams(2:25) \in [-2,2]
% set the first parameter to 1 if gui on
% Dependency: WGCCM_three_link_walker_example/
try
ret = NaN;
cnt = 0;
while isnan(ret) && cnt < 2
ret = walker_main(0,hyperparams(1),'demo',reshape(hyperparams(2:end), [3,8]));
close all;
cnt = cnt+1;
end
if isnan(ret)
    ret = 0;
end

catch MI
    disp(MI)
    allerror = [];
    if exist('error.mat','file') == 2
        load('error.mat');
    end
    allerror = [allerror; hyperparams];
    save('error.mat','allerror','MI')
    ret = 0;
end
