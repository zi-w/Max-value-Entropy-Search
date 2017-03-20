% Copyright (c) 2017 Zi Wang
% Copyright (c) 2014, J.M. Hernandez-Lobato, M.W. Hoffman, Z. Ghahramani
% This function is partially adapted from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
function [ optimum, fval] = globalMaximization(target, xmin, xmax, guesses, useGradient)
% Approximately globally maximize the function target.
% target is a function handle returning the function value (and gradient
% if useGradient is true).
% target can take input of dimension nSample*d if only outputing the values.
%try
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
    % In some situations, the optimization can fail. For example, when the
    % guesses contains the global optimum.
    %disp('fmincon in globalMaximization failed to return a value better than the initialization.')
    [optimum, fval] = fminsearch(neg_target, start);
    fval = -fval;
    if fval < maxVal
        %disp('fminsearch in globalMaximization failed to return a value better than the initialization.')
        
        fval = maxVal;
        optimum = start;
    else
        flag = 0;
        for i = 1:length(xmin)
            if optimum(i) < xmin(i) || optimum(i) > xmax(i)
                flag = 1;
                break;
            end
        end
        if flag
            %disp('fminsearch in globalMaximization failed to search within [xmin, xmax].')
            optimum = start;
            fval = maxVal;
        end
    end
end
% catch ME
%     % In rare situations, the target can return NaN for the gradient, which
%     % causes errors in fmincon.
%     % With 0.2 probability, we sample a random point; with 0.8
%     % probability, we return the best point of the random samples Xgrid.
%     disp(ME.identifier)
%     if rand < 0.8
%         optimum = start;
%         fval = maxVal;
%     else
%         optimum = Xgrid(1,:);
%         fval = y(1);
%     end
% end
end
