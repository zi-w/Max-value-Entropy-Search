clear all; clc;
% add necessary paths
deploy;

% Define function
dx = 10;
xmin = zeros(dx,1);
xmax = ones(dx,1);

f = sample_addGP(dx, dx, xmin, xmax);
options.savefilenm = [];

% Start BO

gpopt(f, xmin, xmax, 200, [], [], options);

% Start BO with add-GP
% add_gpopt(f, xmin, xmax, 200, [], [], options)