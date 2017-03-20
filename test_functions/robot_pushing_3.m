function f = robot_pushing_3(rpos, simu_steps, gpos)
% Returns the distance to goal of an object pushed by a pusher.
% The goal is to minimize f.
% You can define the function to be maximized as following
% gpos = 10 .* rand(1, 2) - 5;
% f = @(x) 5 - robot_pushing_3(x(1:2), x(3), gpos);
% tuning range: 
% xmin = [-5; -5; 1];
% xmax = [5; 5; 30];
f = python('generate_simudata3.py',num2str(rpos(1)),num2str(rpos(2)),num2str(simu_steps),num2str(gpos(1)),num2str(gpos(2)));
f = str2double(f);