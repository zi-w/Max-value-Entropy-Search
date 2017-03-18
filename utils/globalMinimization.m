function [ optimum, fval] = globalMinimization(target, xmin, xmax, guesses)

d = size(xmin, 1);

gridSize = 1000;

Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, 1) .* rand(gridSize, d);
Xgrid = Xgrid';
Xgrid = [ Xgrid guesses ];

y = ones(1, size(Xgrid, 2))*realmax;
for i = 1:size(Xgrid, 2)
    y(i) = target(Xgrid(:,i));
end


[ minValue, minIndex ] = min(y);

start = Xgrid(:,minIndex);

% We optimize starting at the best location

[ optimum, fval] = fmincon(target, start, [], [], [], [], xmin, xmax, [], ...
    optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', 'GradObj', 'off'));
if fval > minValue
    optimum = start;
    fval = minValue;
    disp('failed global opt seq')
end

end
