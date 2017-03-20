function f = robot_pushing_4(rpos, angle, simu_steps, gpos)
% Returns the distance to goal of an object pushed by a pusher.
% The goal is to minimize f.
% You can define the function to be maximized as following
% gpos = 10 .* rand(1, 2) - 5;
% f = @(x) 5 - robot_pushing_4(x(1:2), x(4), x(3), gpos);
% tuning range of x:
% xmin = [-5; -5; 1; 0];
% xmax = [5; 5; 30; 2*pi];
f = python('generate_simudata4.py',num2str(rpos(1)),num2str(rpos(2)),num2str(simu_steps),num2str(gpos(1)),num2str(gpos(2)),num2str(angle));
f = str2double(f);