% Copyright (c) 2017 Zi Wang
function [decomp, hyp] = sampleStructPriors(xx, yy, fixhyp)
% This function optimizes the hyper parameters and the decomposition of
% input dimensions of add-GP.
disp('Start sample struct priors')
dx = size(xx,2);
n_partition = dx;

% if you don't specify fixhyp.l, fixhyp.sigma, or fixhyp.sigma0, it learns
% all the parameters and partition jointly
if isfield(fixhyp, 'l') && isfield(fixhyp, 'sigma') ...
        && isfield(fixhyp, 'sigma0')
    hyp.l = fixhyp.l;
    hyp.sigma = fixhyp.sigma;
    hyp.sigma0 = fixhyp.sigma0;
    [decomp] = learn_partition(xx, yy, hyp, fixhyp, n_partition);
else
    % learn l and sigma and sigma0
    prob_partition = ones(1,n_partition)/n_partition;
    if isfield(fixhyp, 'z')
        decomp = fixhyp.z;
    else
        decomp = sample_categorical(prob_partition, dx);
    end
    % you can increase numIter
    numIter = 2;
    guess_params = [];
    for i = 1:numIter
        assert(max(decomp) <= n_partition)
        nll_func = @(params) compute_nlz_wrap(xx, yy, params, n_partition, decomp);
        bounds = [repmat([0,10], [n_partition, 1]);
            repmat([-5,2], [n_partition, 1]);
            repmat([-10,1], [n_partition, 1])];
        % assume l is the same for all dimensions in each partition
        [best_params, nlltmp] = globalMinimization(nll_func, bounds(:,1), bounds(:,2), guess_params);
        disp(['finished optimize hyp nll=' num2str(nlltmp)])
        guess_params = [guess_params best_params];
        best_params = best_params';
        % params(1:n_partition) defines l in log scale
        l = exp(best_params(1, decomp));
        % params(n_partition+1:2*n_partition) defines sigma in log scale
        sigma = exp(best_params(1, n_partition+1:2*n_partition));
        % params(2*n_partition+1:3*n_partition) defines sigma0 in log scale
        sigma0 = exp(best_params(1, 2*n_partition+1:3*n_partition));
        hyp.l = l;
        hyp.sigma = sigma;
        hyp.sigma0 = sigma0;
        % learn z
        [decomp] = learn_partition(xx, yy, hyp, fixhyp, n_partition);
        fixhyp.z = decomp;
    end
    
end

end

function [decomp] = learn_partition(xx, yy, hyp, fixhyp, n_partition)
if isfield(fixhyp, 'decomp')
    decomp = fixhyp.decomp;
    return;
end
%n_partition = fixhyp.n_partition;
% Number of iterations for gibbs sampling.
N_gibbs = 10;
gibbs_iter = floor(N_gibbs/2);
% Limit the number of dimensions in each group because BO cannot do well in
% high dimensions.
dim_limit = 5;

% You may want to set maxNdata to be small to speed up this function.
maxNdata = 750;
dx = size(xx,2);
Nidx = min(maxNdata, size(xx,1));
xx = xx(1:Nidx,:);
yy = yy(1:Nidx);

% This dirichlet hyper parameter is something one can tune.
hyp_dirichlet = ones(1, n_partition) *1;
prob_partition = hyp_dirichlet ./ sum(hyp_dirichlet);
if isfield(fixhyp, 'z')
    z = fixhyp.z;
    z_best = z;
else
    z = sample_categorical(prob_partition, dx);
    
    z_best = 0;
end

minnlz = realmax;

% best z without constraints
minnlz_nocon = realmax;
%% start gibbs iteration
for i = 1:N_gibbs
    for d = 1:dx
        
        if i < gibbs_iter % in the first two iterations, no limit
            log_prob = zeros(1, n_partition);
            nlz = zeros(1, n_partition);
            all_cat = unique(z);
            other_cat = 1:n_partition;
            % calculate the categorical distribution
            for a = all_cat
                other_cat = other_cat(other_cat ~= a);
                z(d) = a;
                nlz(a) =  compute_nlz(xx, yy, hyp, z);
                
                log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
            end
            if ~isempty(other_cat)
                z(d) = other_cat(1);
                nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
            end
            
            [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
            
            if minnlz_nocon > nlz(z(d))
                z_best = z;
                minnlz_nocon = nlz(z(d));
            end
            if i == gibbs_iter-1
                z = z_best;
            end
        else % limit #partitions
            log_prob = zeros(1, n_partition);
            nlz = zeros(1, n_partition);
            all_cat = unique(z);
            other_cat = 1:n_partition;
            for a = all_cat
                z(d) = -1;
                other_cat = other_cat(other_cat ~= a);
                if sum(z == a) < dim_limit
                    
                    z(d) = a;
                    nlz(a) =  compute_nlz(xx, yy, hyp, z);
                    
                    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                else
                    log_prob(a) = -realmax;
                    nlz(a) = realmax;
                end
            end
            if ~isempty(other_cat)
                z(d) = other_cat(1);
                nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                
                log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
            end
            
            [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
            
            if minnlz > nlz(z(d))  && i > gibbs_iter
                z_best = z;
                minnlz = nlz(z(d));
            end
            
        end
    end
    
end
decomp = z_best;

end