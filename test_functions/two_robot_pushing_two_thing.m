function f = two_robot_pushing_two_thing(rpos, rvel, simu_steps, angle, rpos2, ...
    rvel2, simu_steps2, angle2, rtor, rtor2, gpos, gpos2)
% Returns the sum of distances to goal of two objects pushed by two pushers.
% The goal is to minimize f.
% You can define the function to be maximized as following
% gpos1 = [4 3.5];
% gpos2 = [-4 -3.5];
% f = @(x) 10 - two_robot_pushing_two_thing(x(1:2), x(3:4), x(5), x(6), ...
%    x(7:8), x(9:10), x(11), x(12), x(13), x(14), gpos1, gpos2);
% tuning range of x:
% xmin = [-5; -5; -10; -10; 2; 0; -5; -5; -10; -10; 2; 0; -5; -5];
% xmax = [5; 5; 10; 10; 30; 2*pi; 5;5;10;10;30;2*pi; 5; 5];
f= nan;
cnt = 0;
while isnan(f) && cnt < 10
    f = python('generate_simudata_2robot2thing.py',num2str(rpos(1)),num2str(rpos(2)), ...
        num2str(rvel(1)),num2str(rvel(2)),num2str(simu_steps),num2str(angle), ...
        num2str(rpos2(1)),num2str(rpos2(2)), ...
        num2str(rvel2(1)),num2str(rvel2(2)),num2str(simu_steps2),num2str(angle2), ...
        num2str(rtor), num2str(rtor2), ...
        num2str(gpos(1)),num2str(gpos(2)), num2str(gpos2(1)),num2str(gpos2(2)));
    f = str2double(f);
    cnt = cnt + 1;
end
if isnan(f)
    f = 0;
end