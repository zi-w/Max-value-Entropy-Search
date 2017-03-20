% Copyright (c) 2017 Zi Wang
function Xgrid = rand_sample_interval(xmin, xmax, gridSize)
    % Uniformly randomly sample gridSize number of points in intervals
    % defined by xmin and xmax.
    % xmin is a d*1 vector specifying the lower bound;
    % xmax is a d*1 vector specifying the upper bound;
    % gridSize is the number of samples to return;
    % Xgrid is a gridSize*d matrix.
    
	d = size(xmin, 1);
    
	Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, 1) .* rand(gridSize, d);