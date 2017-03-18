% Author: Zi Wang
% See also: gpopt.m
function results = add_gpopt(objective, xmin, xmax, T, initx, inity, options)
% This function maximizes the function objective using BO with add-GP and 
% returns results as a cell of size 7, including the inferred argmax points 
% (guesses),the function values of the inferred argmax points (guessvals), the
% evaluated points (xx), the function values of the evaluated points
% (yy), the runtime to choose the points (choose_time) and extra time of
% inferring the argmax (extra_time).
% objective is a function handle;
% xmin is a column vector indicating the lower bound of the search space;
% xmax is a column vector indicating the upper bound of the search space;
% T is the number of sequential evaluations of the function;
% initx, inity are the initialization of observed inputs and their values.

if nargin <= 6
    options = struct();
end
if ~isfield(options, 'restart') options.restart = 0; end
if ~isfield(options, 'bo_method'); options.bo_method = 'mes-g'; end
if ~isfield(options, 'savefilenm'); options.savefilenm = []; end
% When testing synthetic functions, one can add noise to the output.
if ~isfield(options, 'noiselevel'); options.noiselevel = 0; end
if ~isfield(options, 'nK'); options.nK = 1; end
if ~isfield(options, 'nFeatures'); options.nFeatures = 10000; end
if ~isfield(options, 'seed'); options.seed = 42; end
% Set random seed
s = RandStream('mcg16807','Seed', options.seed);
if ~isfield(options, 'learn_interval'); options.learn_interval = 10; end
RandStream.setGlobalStream(s);

% Set options.restart = 1 to use the saved results and run more iterations
if options.restart && exist(savefilenm,'file') ~= 2
    options.restart = 0;
end
if options.restart == 0
    if isempty(initx)
        % initialize xx,yy with at least one pair of intx, inty
        initx = rand_sample_interval(xmin, xmax, 1);
        inity = objective(initx);
    end
    xx = initx;
    yy = inity;
    guesses = initx;
    guessvals = inity;
    choose_time = []; % elapsed time to choose where to evaluate
    extra_time = []; % elapsed time to optimize mean function, hyper-parameters
    tstart = 0;
else
    restart_file = load(options.savefilenm);
    guesses      =  restart_file.results{1};
    guessvals    =  restart_file.results{2};
    xx           =  restart_file.results{3};
    yy           =  restart_file.results{4};
    choose_time  =  restart_file.results{5};
    extra_time   =  restart_file.results{6};
    t            =  restart_file.results{7};
    z            =  restart_file.results{8};
    tstart = t;
    if tstart >= T
        return
    end
end


%% start optimization
for t = tstart+1 : T
    
    tic
    if mod(t, fixhyp.learn_interval) == 1
        % learn structure
        [z, hyp, minnlz ] = sampleStructPriors(xx, yy, fixhyp);
        fixhyp.z = z;
    end
    extra_time = [extra_time; toc];
    
    tic
    % Calculate and inverse the gram matrix
    KernelMatrixInv = cell(1);
    KernelMatrix = compute_gram(xx, hyp, 1, z);
    KernelMatrixInv{1} = chol2invchol(KernelMatrix);

    all_cat = unique(z);
    if strcmp(bo_method, 'MES-R')
        % compute the max-values maxes of size nZ x nK
        subnFeatures = ceil(optioins.nFeatures/length(all_cat));
        maxes = add_sampleMaximumValues(1, nK, xx, yy, hyp.sigma0, hyp.sigma, hyp.l, xmin, xmax, subnFeatures, z, all_cat);
    end
    xnext = zeros(1, size(xx,2));
    % Start optimization group by group
    for i = 1:length(all_cat)
        coords = (z==all_cat(i));
        xx_sub = xx(:,coords);
        xmin_sub = xmin(coords);
        xmax_sub = xmax(coords);
        l = hyp.l(:,coords);
        sigma = hyp.sigma(:,all_cat(i));
        sigma0 = hyp.sigma0(:,all_cat(i));
        if fixhyp.usediscrete
            fixhyp.cur_discrete = fixhyp.discrete(coords,:);
            fixhyp.cur_grid = fixhyp.xgrid{i};
        end
        if strcmp(bo_method, 'UCB')
            % The parameter is set to be the same as add-GP-UCB 
            % (Kandasamy et al. ICML 2015)
            alpha = 1;
            beta = sqrt(size(xx_sub,2)*log(2*t)/5);
            optimum = ucb_choose(xx_sub, yy, KernelMatrixInv, [], ...
                sigma0, sigma, l, xmin_sub, xmax_sub, alpha, beta);
   
        elseif strcmp(bo_method, 'MES-R')
            maxes_sub = maxes(i,:);
            optimum = add_mesr_choose(maxes_sub, xx_sub, yy, ...
                KernelMatrixInv,  sigma0, sigma, l, xmin_sub, xmax_sub);

        elseif strcmp(bo_method, 'MES-G')
            optimum = add_mesg_choose(nK, xx_sub, yy, KernelMatrixInv, sigma0, sigma, l, xmin_sub, xmax_sub, theta_shift, islog, logfileID);

        elseif strcmp(bo_method, 'EST')
            optimum = add_est_choose(nM, xx_sub, yy, KernelMatrixInv, sigma0, sigma, l, xmin_sub, xmax_sub);
        end
        xnext(coords) = optimum;
    end
    choose_time = [choose_time; toc];
    
    xx = [ xx ; xnext ];
    yy = [ yy ; objective(xnext) + rand(1) * noiselevel];
    
    disp([num2str(t) ': ' 'val=' num2str(yy(end,:))])
    % save result every a few iterations
    if ~isempty(options.savefilenm) && mod(t,10) == 0
        results{1} = guesses;
        results{2} = guessvals;
        results{3} = xx;
        results{4} = yy;
        results{5} = choose_time;
        results{6} = extra_time;
        results{7} = t;
        results{8} = z;
        save(savefilenm, 'results');
    end
end