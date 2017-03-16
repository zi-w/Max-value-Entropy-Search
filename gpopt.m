% Author: Zi Wang
% This function is adapted from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
% https://bitbucket.org/jmh233/codepesnips2014
function results = gpopt(objective, xmin, xmax, T, initx, inity, options)
% xmin xmax are column vectors
if nargin <= 6
    options = struct();
end
if ~isfield(options, 'restart') options.restart = 0; end
if ~isfield(options, 'bo_method'); options.bo_method = 'mes-g'; end
if ~isfield(options, 'savefilenm'); options.savefilenm = []; end
if ~isfield(options, 'noiselevel'); options.noiselevel = 0; end
if ~isfield(options, 'nM'); options.nM = 5; else nM = 1; end
if ~isfield(options, 'nK'); options.nK = 1; end
if ~isfield(options, 'nFeatures'); options.nFeatures = 1000; end
if ~isfield(options, 'seed'); options.seed = 1000; end
% Set random seed
s = RandStream('mcg16807','Seed', options.seed);
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
    tstart = t;
    if tstart >= T
        return
    end
end

% We sample from the posterior distribution of the hyper-parameters

[ l, sigma, sigma0 ] = sampleHypers(xx, yy, nM, fixhyp);
KernelMatrixInv = cell(1, nM);
for j = 1 : nM
    KernelMatrix = computeKmm(xx, l(j,:)', sigma(j), sigma0(j));
    KernelMatrixInv{ j } = chol2invchol(KernelMatrix);
end

results = cell(1,7);
for t = tstart+1 : T
    
    tic
    
    if strcmp(options.bo_method, 'MES-R')
        optimum = mesr_choose(nM, nK, xx, yy, KernelMatrixInv, ...
            guesses, sigma0, sigma, l, xmin, xmax, nFeatures);
    elseif strcmp(options.bo_method, 'MES-G')
        optimum = mesg_choose(nM, nK, xx, yy, KernelMatrixInv, ...
            guesses, sigma0, sigma, l, xmin, xmax);
    elseif strcmp(options.bo_method, 'EI')
        optimum = ei_choose(xx, yy, KernelMatrixInv, guesses, ...
            sigma0, sigma, l, xmin, xmax);
    elseif strcmp(options.bo_method, 'PI')
        optimum = pi_choose(nM, nK, xx, yy, KernelMatrixInv, guesses, ...
            sigma0, sigma, l, xmin, xmax, nFeatures, m0);
    elseif strcmp(options.bo_method, 'UCB')
        alpha = 1;
        beta = (2*log(t^2*2*pi^2/(3*0.01)) + 2*length(xmin)*log(t^2*...
            length(xmin)*max(xmax-xmin)*(log(4*length(xmin)/0.01))^0.5))^0.5;
        optimum = ucb_choose(xx, yy, KernelMatrixInv, guesses, ...
            sigma0, sigma, l, xmin, xmax, alpha, beta);
    elseif strcmp(options.bo_method, 'est')
        optimum = est_choose(nM, xx, yy, KernelMatrixInv, guesses, ...
            sigma0, sigma, l, xmin, xmax);
    end
    
    
    xx = [ xx ; optimum ];
    
    yy = [ yy ; objective(optimum)+ randn(1)*noiselevel];
    
    if mod(t, fixhyp.learn_interval) == 0
        % We sample from the posterior distribution of the hyper-parameters
        [ l, sigma, sigma0 ] = sampleHypers(xx, yy, nM, fixhyp);
    end
    % We update the kernel matrix on the samples
    
    KernelMatrixInv = cell(1, nM);
    for j = 1 : nM
        KernelMatrix = computeKmm(xx, l(j,:)', sigma(j), sigma0(j));
        KernelMatrixInv{ j } = chol2invchol(KernelMatrix);
    end
    
    choose_time = [choose_time; toc];
    
    
    tic
    % We optimize the posterior mean of the GP
    f = @(x) posteriorMean(x, xx, yy, KernelMatrixInv, l, sigma);
    
    optimum = globalMaximization(f, xmin, xmax, guesses);
    extra_time = [extra_time; toc];
    
    guesses = [ guesses ; optimum ];
    guessvals = [guessvals; objective(optimum)];
    
    disp(['elapsed choosing time is ' num2str(choose_time(end))])
    disp([num2str(t) ': tested ' num2str(xx(end,:)) '; val=' num2str(yy(end,:)) ...
        '; guess ' num2str(optimum) ' val=' num2str(guessvals(end))])
    
    % Save result every a few iterations
    if  ~isempty(options.savefilenm) && mod(t, 10) == 0
        results{1} = guesses;
        results{2} = guessvals;
        results{3} = xx;
        results{4} = yy;
        results{5} = choose_time;
        results{6} = extra_time;
        results{7} = t;
        save(savefilenm, 'results');
    end
end