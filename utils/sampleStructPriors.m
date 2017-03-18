function [decomp, hyp, minnlz] = sampleStructPriors(xx, yy, nM, fixhyp)
disp('Start sample struct priors')
dx = size(xx,2);
n_partition = fixhyp.n_partition;
% todo: learn and sample l, sigma, sigma0

% if you don't specify fixhyp.l, fixhyp.sigma, or fixhyp.sigma0, it learns
% all the parameters and partition jointly
if isfield(fixhyp, 'l') && isfield(fixhyp, 'sigma') && isfield(fixhyp, 'sigma0')
    hyp.l = fixhyp.l;
    hyp.sigma = fixhyp.sigma;
    hyp.sigma0 = fixhyp.sigma0;
    [decomp, minnlz] = learn_partition(xx, yy, hyp, fixhyp);
else
    % learn l and sigma and sigma0
    prob_partition = ones(1,n_partition)/n_partition;%[0.5 0.5];
    if isfield(fixhyp, 'z')
        decomp = fixhyp.z;
    else
        decomp = sample_categorical(prob_partition, dx);
    end
    numIter = 1;
    guess_params = [];
    for i = 1:numIter
        disp(['sample struct prior iter ' num2str(i)])
        % learn l, sigma, sigma0 via direct
        assert(max(decomp) <= n_partition)
        nll_func = @(params) compute_nlz_wrap(xx, yy, params, n_partition, decomp, fixhyp.flag_discrete);
        % direct is extremely slow
        %Problem.f = nll_func;
        %opts.maxits = 10;
        bounds = [repmat([0,10], [n_partition, 1]);
            repmat([-5,2], [n_partition, 1]);
            repmat([-10,1], [n_partition, 1])];
        %[~, best_params] = diRect(Problem, bounds, opts);
        % params(1:n_partition) defines l in log scale
        % assume l is the same for all dimensions in each partition
        [best_params, nlltmp] = globalOptimization_seq(nll_func, bounds(:,1), bounds(:,2), guess_params);
        disp(['finished optimize hyp nll=' num2str(nlltmp)])
        guess_params = [guess_params best_params];
        best_params = best_params';
        l = exp(best_params(1, decomp));
        l(fixhyp.flag_discrete > 0) = -l(fixhyp.flag_discrete > 0);
        %if any(l > 0)
        %    disp('in sample struct')
        %    keyboard;
        %end
        % params(n_partition+1:2*n_partition) defines sigma in log scale
        sigma = exp(best_params(1, n_partition+1:2*n_partition));
        % params(2*n_partition+1:3*n_partition) defines sigma0 in log scale
        sigma0 = exp(best_params(1, 2*n_partition+1:3*n_partition));
        hyp.l = l;
        hyp.sigma = sigma;
        hyp.sigma0 = sigma0;
        % learn z
        [decomp, minnlz] = learn_partition(xx, yy, hyp, fixhyp);
        disp(['finished learn partition nll=' num2str(minnlz)])
        fixhyp.z = decomp;
    end
    
end

end

function [decomp, minnlz] = learn_partition(xx, yy, hyp, fixhyp)
n_partition = fixhyp.n_partition;
N_gibbs = 8;
maxNdata = 750;
dx = size(xx,2);
if fixhyp.decomp_learn == 0 || fixhyp.decomp_learn == 3 || fixhyp.decomp_learn == 5
    decomp = fixhyp.decomp;
    minnlz = 0;
elseif fixhyp.decomp_learn == 1 || fixhyp.decomp_learn == 6% use gibbs sampling to learn decomp
    % choose at most 350 points
    %
    %idxs = randperm(size(xx,1));
    %xx = xx(idxs(1:Nidx),:);
    %yy = yy(idxs(1:Nidx));
    Nidx = min(maxNdata, size(xx,1));
    xx = xx(1:Nidx,:);
    yy = yy(1:Nidx);
    hyp_dirichlet = ones(1, n_partition) ;%* 0.01;
    prob_partition = hyp_dirichlet ./ sum(hyp_dirichlet);
    if isfield(fixhyp, 'z')
        z = fixhyp.z;
    else
        z = sample_categorical(prob_partition, dx);
    end
    
    %nlzs = zeros(dx, N_gibbs);
    minnlz = compute_nlz(xx, yy, hyp, z);
    curnlz = minnlz;
    z_best = z;
    
    %% start gibbs iteration
    for i = 1:N_gibbs
        for d = 1:dx
            % calculate the categorical distribution
            log_prob = zeros(1, n_partition);
            nlz = zeros(1, n_partition);
            all_cat = unique(z);
            orizd = z(d);
            other_cat = 1:n_partition;
            
            for a = all_cat%1:n_partition
                other_cat = other_cat(other_cat ~= a);
                if a == orizd
                    z(d) = a;
                    nlz(a) = curnlz;
                    %assert(abs(curnlz - compute_nlz(xx, yy, hyp, z))<0.001)
                else
                    z(d) = a;
                    nlz(a) =  compute_nlz(xx, yy, hyp, z);
                end
                log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                %- 0.5*yy'/(K + hyp.sigma0(a)*eye(n_x))*yy - 0.5*log(det(K + hyp.sigma0(a)*eye(n_x)));
            end
            if ~isempty(other_cat)
                z(d) = other_cat(1);
                nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                %for a = other_cat
                %    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                %end
                log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
            end
            %prob = -min(log_prob) + log_prob;
            % sample from this distribution
            %prob = exp(prob);
            %prob = prob/sum(prob);
            %if rand < 0
            %    [~,z(d)] = max(prob);
            %else
            %    z(d) = sample_categorical(prob, 1);
            %end
            [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
            %nlzs(d,i) = nlz(z(d));
            %disp(prob(z(d)))
            curnlz = nlz(z(d));
            if minnlz > nlz(z(d))
                z_best = z;
                minnlz = nlz(z(d));%nlzs(d,i);
            end
        end
        
    end
    decomp = z_best;
elseif fixhyp.decomp_learn == 7 % limit number of dimensions in each group
    
    dim_limit = 2;
    gibbs_iter = N_gibbs-4; % number of iterations running full gibbs
    % choose at most 350 points
    %
    %idxs = randperm(size(xx,1));
    %xx = xx(idxs(1:Nidx),:);
    %yy = yy(idxs(1:Nidx));
    Nidx = min(maxNdata, size(xx,1));
    xx = xx(1:Nidx,:);
    yy = yy(1:Nidx);
    hyp_dirichlet = ones(1, n_partition) *1;%* 0.01;
    prob_partition = hyp_dirichlet ./ sum(hyp_dirichlet);
    if isfield(fixhyp, 'z')
        z = fixhyp.z;
        minnlz = compute_nlz(xx, yy, hyp, z);
        curnlz = minnlz;
        z_best = z;
    else
        z = sample_categorical(prob_partition, dx);
        minnlz = compute_nlz(xx, yy, hyp, z);
        curnlz = minnlz;
        
        z_best = 0;
    end
    
    %nlzs = zeros(dx, N_gibbs);
    minnlz = realmax;
    
    % best z without constraints
    z_best_nocon = 0;
    minnlz_nocon = realmax;
    %% start gibbs iteration
    for i = 1:N_gibbs
        for d = 1:dx
            % calculate the categorical distribution
            if i < gibbs_iter % in the first two iterations, no limit
                log_prob = zeros(1, n_partition);
                nlz = zeros(1, n_partition);
                all_cat = unique(z);
                orizd = z(d);
                other_cat = 1:n_partition;
                for a = all_cat%1:n_partition
                    other_cat = other_cat(other_cat ~= a);
                    %if a == orizd
                    %    z(d) = a;
                    %    nlz(a) = curnlz;
                    %    assert(abs(curnlz - compute_nlz(xx, yy, hyp, z))<0.001)
                    %else
                        z(d) = a;
                        nlz(a) =  compute_nlz(xx, yy, hyp, z);
                    %end
                    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %- 0.5*yy'/(K + hyp.sigma0(a)*eye(n_x))*yy - 0.5*log(det(K + hyp.sigma0(a)*eye(n_x)));
                end
                if ~isempty(other_cat)
                    z(d) = other_cat(1);
                    nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                    %for a = other_cat
                    %    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %end
                    log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
                end
                %prob = -min(log_prob) + log_prob;
                % sample from this distribution
                %prob = exp(prob);
                %prob = prob/sum(prob);
                %if rand < 0
                %    [~,z(d)] = max(prob);
                %else
                %    z(d) = sample_categorical(prob, 1);
                %end
                [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
                %nlzs(d,i) = nlz(z(d));
                %disp(prob(z(d)))
                
                curnlz = nlz(z(d));
                if minnlz_nocon > nlz(z(d))
                    z_best_nocon = z;
                    minnlz_nocon = nlz(z(d));
                end
                if i == gibbs_iter-1
                    z = z_best_nocon;
                    curnlz = minnlz_nocon;
                end
            else % limit #partitions
                log_prob = zeros(1, n_partition);
                nlz = zeros(1, n_partition);
                all_cat = unique(z);
                orizd = z(d);
                other_cat = 1:n_partition;
                for a = all_cat%1:n_partition
                    z(d) = -1;
                    other_cat = other_cat(other_cat ~= a);
                    if sum(z == a) < dim_limit
                        %if a == orizd
                        %    z(d) = a;
                        %    nlz(a) = curnlz;
                        %    assert(abs(curnlz - compute_nlz(xx, yy, hyp, z))<0.001)
                        %else
                            z(d) = a;
                            nlz(a) =  compute_nlz(xx, yy, hyp, z);
                        %end
                        log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    else
                        log_prob(a) = -realmax;
                        nlz(a) = realmax;
                    end
                    %- 0.5*yy'/(K + hyp.sigma0(a)*eye(n_x))*yy - 0.5*log(det(K + hyp.sigma0(a)*eye(n_x)));
                end
                if ~isempty(other_cat)
                    z(d) = other_cat(1);
                    nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                    %for a = other_cat
                    %    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %end
                    log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
                end
                %prob = -min(log_prob) + log_prob;
                % sample from this distribution
                %prob = exp(prob);
                %prob = prob/sum(prob);
                %if rand < 0
                %    [~,z(d)] = max(prob);
                %else
                %    z(d) = sample_categorical(prob, 1);
                %end
                [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
                %nlzs(d,i) = nlz(z(d));
                %disp(prob(z(d)))
                curnlz = nlz(z(d));
                if minnlz > nlz(z(d))  && i > gibbs_iter
                    z_best = z;
                    minnlz = nlz(z(d));%nlzs(d,i);
                end
                
            end
        end
        
    end
    decomp = z_best;
elseif fixhyp.decomp_learn == 8 % limit number of dimensions in each group
    
    dim_limit = 3;
    gibbs_iter = N_gibbs-4; % number of iterations running full gibbs
    % choose at most 350 points
    %
    %idxs = randperm(size(xx,1));
    %xx = xx(idxs(1:Nidx),:);
    %yy = yy(idxs(1:Nidx));
    Nidx = min(maxNdata, size(xx,1));
    xx = xx(1:Nidx,:);
    yy = yy(1:Nidx);
    hyp_dirichlet = ones(1, n_partition) *1;%* 0.01;
    prob_partition = hyp_dirichlet ./ sum(hyp_dirichlet);
    if isfield(fixhyp, 'z')
        z = fixhyp.z;
        minnlz = compute_nlz(xx, yy, hyp, z);
        curnlz = minnlz;
        z_best = z;
    else
        z = sample_categorical(prob_partition, dx);
        minnlz = compute_nlz(xx, yy, hyp, z);
        curnlz = minnlz;
        
        z_best = 0;
    end
    
    %nlzs = zeros(dx, N_gibbs);
    minnlz = realmax;
    
    % best z without constraints
    z_best_nocon = 0;
    minnlz_nocon = realmax;
    %% start gibbs iteration
    for i = 1:N_gibbs
        for d = 1:dx
            % calculate the categorical distribution
            if i < gibbs_iter % in the first two iterations, no limit
                log_prob = zeros(1, n_partition);
                nlz = zeros(1, n_partition);
                all_cat = unique(z);
                orizd = z(d);
                other_cat = 1:n_partition;
                for a = all_cat%1:n_partition
                    other_cat = other_cat(other_cat ~= a);
                    %if a == orizd
                    %    z(d) = a;
                    %    nlz(a) = curnlz;
                    %    assert(abs(curnlz - compute_nlz(xx, yy, hyp, z))<0.001)
                    %else
                        z(d) = a;
                        nlz(a) =  compute_nlz(xx, yy, hyp, z);
                    %end
                    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %- 0.5*yy'/(K + hyp.sigma0(a)*eye(n_x))*yy - 0.5*log(det(K + hyp.sigma0(a)*eye(n_x)));
                end
                if ~isempty(other_cat)
                    z(d) = other_cat(1);
                    nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                    %for a = other_cat
                    %    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %end
                    log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
                end
                %prob = -min(log_prob) + log_prob;
                % sample from this distribution
                %prob = exp(prob);
                %prob = prob/sum(prob);
                %if rand < 0
                %    [~,z(d)] = max(prob);
                %else
                %    z(d) = sample_categorical(prob, 1);
                %end
                [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
                %nlzs(d,i) = nlz(z(d));
                %disp(prob(z(d)))
                
                curnlz = nlz(z(d));
                if minnlz_nocon > nlz(z(d))
                    z_best_nocon = z;
                    minnlz_nocon = nlz(z(d));
                end
                if i == gibbs_iter-1
                    z = z_best_nocon;
                    curnlz = minnlz_nocon;
                end
            else % limit #partitions
                log_prob = zeros(1, n_partition);
                nlz = zeros(1, n_partition);
                all_cat = unique(z);
                orizd = z(d);
                other_cat = 1:n_partition;
                for a = all_cat%1:n_partition
                    z(d) = -1;
                    other_cat = other_cat(other_cat ~= a);
                    if sum(z == a) < dim_limit
                        %if a == orizd
                        %    z(d) = a;
                        %    nlz(a) = curnlz;
                        %    assert(abs(curnlz - compute_nlz(xx, yy, hyp, z))<0.001)
                        %else
                            z(d) = a;
                            nlz(a) =  compute_nlz(xx, yy, hyp, z);
                        %end
                        log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    else
                        log_prob(a) = -realmax;
                        nlz(a) = realmax;
                    end
                    %- 0.5*yy'/(K + hyp.sigma0(a)*eye(n_x))*yy - 0.5*log(det(K + hyp.sigma0(a)*eye(n_x)));
                end
                if ~isempty(other_cat)
                    z(d) = other_cat(1);
                    nlz(other_cat) = compute_nlz(xx, yy, hyp, z);
                    %for a = other_cat
                    %    log_prob(a) = log(sum(z==a) + hyp_dirichlet(a)) - nlz(a);
                    %end
                    log_prob(other_cat) = log(1 + hyp_dirichlet(other_cat)) - nlz(other_cat);
                end
                %prob = -min(log_prob) + log_prob;
                % sample from this distribution
                %prob = exp(prob);
                %prob = prob/sum(prob);
                %if rand < 0
                %    [~,z(d)] = max(prob);
                %else
                %    z(d) = sample_categorical(prob, 1);
                %end
                [~,z(d)] = max(log_prob - log(-log(rand(1, n_partition))));
                %nlzs(d,i) = nlz(z(d));
                %disp(prob(z(d)))
                curnlz = nlz(z(d));
                if minnlz > nlz(z(d))  && i > gibbs_iter
                    z_best = z;
                    minnlz = nlz(z(d));%nlzs(d,i);
                end
                
            end
        end
        
    end
    decomp = z_best;
elseif fixhyp.decomp_learn == 2 % partial optimize by random sampling
    %Nidx = min(maxNdata, size(xx,1));
    %idxs = randperm(size(xx,1));
    %xx = xx(idxs(1:Nidx),:);
    %yy = yy(idxs(1:Nidx));
    Nidx = min(maxNdata, size(xx,1));
    xx = xx(1:Nidx,:);
    yy = yy(1:Nidx);
    prob_partition = ones(1,n_partition)/n_partition;%[0.5 0.5];
    if isfield(fixhyp, 'z')
        z = fixhyp.z;
    else
        z = sample_categorical(prob_partition, dx);
    end
    minnlz = compute_nlz(xx, yy, hyp, z);
    z_rand_best = z;
    for i = 1:dx*N_gibbs%*n_partition
        z = sample_categorical(prob_partition, dx);
        nlzs2 = compute_nlz(xx, yy, hyp, z);
        if nlzs2 < minnlz
            z_rand_best = z;
            minnlz = nlzs2;
        end
    end
    decomp = z_rand_best;
elseif fixhyp.decomp_learn == 4 % partial optimize by random sampling
    % set # of trials to be the same as in ADD GP UCB
    numTrials = 5;
    %Nidx = min(maxNdata, size(xx,1));
    %idxs = randperm(size(xx,1));
    %xx = xx(idxs(1:Nidx),:);
    %yy = yy(idxs(1:Nidx));
    Nidx = min(maxNdata, size(xx,1));
    xx = xx(1:Nidx,:);
    yy = yy(1:Nidx);
    prob_partition = ones(1,n_partition)/n_partition;%[0.5 0.5];
    if isfield(fixhyp, 'z')
        z = fixhyp.z;
    else
        z = sample_categorical(prob_partition, dx);
    end
    minnlz = compute_nlz(xx, yy, hyp, z);
    z_rand_best = z;
    for i = 1:numTrials%dx*N_gibbs*n_partition
        z = sample_categorical(prob_partition, dx);
        nlzs2 = compute_nlz(xx, yy, hyp, z);
        if nlzs2 < minnlz
            z_rand_best = z;
            minnlz = nlzs2;
        end
    end
    decomp = z_rand_best;
end
end