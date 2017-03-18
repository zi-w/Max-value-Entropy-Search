deploy;

expi = 2;
% Define function
fun = load(['../gorf/test_funcs/3d-' num2str(expi) '.mat']);

options.isfix = 1;
options.l = fun.l;
options.sigma = fun.sigma;
options.sigma0 = fun.sigma0;
options.bo_method = 'MES-R';
options.noiselevel = sqrt(fun.sigma0);
options.nM = 1;
options.nK = 10;
options.savefilenm = ['results3/' options.bo_method '_' num2str(expi) '.mat'];

% Start BO

gpopt(fun.f, fun.xlimit(1,:)', fun.xlimit(2,:)', 200, fun.xx, fun.yy, options);