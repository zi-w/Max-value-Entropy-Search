deploy;

expi = 2;
% Define function
fun = load(['../gorf/test_funcs/3d-' num2str(expi) '.mat']);

options.isfix = 1;
%options.l = fun.l;
%options.sigma = fun.sigma;
%options.sigma0 = fun.sigma0;
options.bo_method = 'add-EST';
options.noiselevel = sqrt(fun.sigma0);
options.nM = 1;
options.nK = 1;
options.savefilenm = ['results3/' options.bo_method '_' num2str(expi) '.mat'];

% Start BO

%gpopt(fun.f, fun.xlimit(1,:)', fun.xlimit(2,:)', 200, fun.xx, fun.yy, options);

% Start BO with add-GP
%options.decomp = [1 1 2];
%options.sigma = repmat(fun.sigma,[1, 2]);
%options.sigma0 = repmat(fun.sigma0, [1, 2]);
add_gpopt(fun.f, fun.xlimit(1,:)', fun.xlimit(2,:)', 200, fun.xx, fun.yy, omaxptions)