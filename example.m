% Copyright (c) 2017 Zi Wang
clear all; clc;
% add necessary paths
deploy;

% Define function
dx = 3;
xmin = zeros(dx,1);
xmax = ones(dx,1);

f = sample_addGP(dx, dx, xmin, xmax);

% Save the file to a path 
options.savefilenm = [];

% Choose BO methods
% options.bo_method = 'Add-MES-G';

% Set the number of maximums to sample
options.nK = 5;

% Set the GP hyper-parameters if you would like to fix them.
% Comment the following 3 lines out if you would like to learn them.
options.l = ones(1,dx)*50;
options.sigma = ones(1, dx)*5;
options.sigma0 = 0.0001*ones(1, dx);

% Start BO
% Set the number of GP hyper-parameter settings to be sampled
% options.nM = 10;
gpopt(f, xmin, xmax, 200, [], [], options);

% Start BO with add-GP
% add_gpopt does not support sampling multiple hyper-parameter settings.
% see sampleStructPriors.m for more details on learning hyper-parameters
% and the additive structure of the function.
% The additive learning strategy is based on the paper
% Wang, Zi and Li, Chengtao and Jegelka, Stefanie and Kohli, Pushmeet. 
% Batched High-dimensional Bayesian Optimization via Structural Kernel 
% Learning. arXiv preprint arXiv:1703.01973.
% add_gpopt(f, xmin, xmax, 200, [], [], options)
