% Author: Zi Wang
% This function is adapted from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
function [ optimum, fval] = globalMaximization(target, xmin, xmax, guesses, useGradient)
% Approximately globally maximize the function target.
% target is a function handle returning the function value (and gradient
% if useGradient is true).
% target can take input of dimension nSample*d if only outputing the values.
try
    % Sample a random set of inputs and find the best one to initialize
    % fmincon for optimization
    if nargin == 3
        guesses = [];
    end
    if nargin == 4
        useGradient = 1;
    end
    gridSize = 10000;
    Xgrid = rand_sample_interval(xmin, xmax, gridSize);
    Xgrid = [ Xgrid ; guesses ];
    
    y = target(Xgrid);
    [maxVal, maxIdx] = max(y);
    
    start = Xgrid(maxIdx,:);
    
    % Wrap the target to its negative values
    gradObj = 'on';
    if useGradient
        neg_target = @(x) negative_wrapper(target, x);
    else
        neg_target = @(x) -target(x);
        gradObj = 'off';
    end
    [ optimum, fval] = fmincon(neg_target, start, [], [], [], [], xmin, xmax, [], ...
        optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', 'GradObj', gradObj));
    fval = -fval;
    
    if fval < maxVal
        % In rare situations, the optimization can fail
        fval = maxVal;
        optimum = start;
    end
catch ME
    % In rare situations, the target can return NaN for the gradient, which
    % causes errors in fmincon.
    % With 0.2 probability, we sample a random point; with 0.8
    % probability, we return the best point of the random samples Xgrid.
    disp(ME.identifier)
    if rand < 0.8
        optimum = start;
        fval = maxVal;
    else
        optimum = Xgrid(1,:);
        fval = y(1);
    end
end
end
