function optimum = add_mesg_choose(nK, xx, yy, KernelMatrixInv, ...
    sigma0, sigma, l, xmin, xmax)
% This function returns the next evaluation point using add-MES-G.
% nK is the number of sampled maximum values.
% xx, yy are the current observations.
% KernelMatrixInv is the gram maxtrix inverse.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.

gridSize = 10000;
d = size(xmin, 1);

Xgrid = repmat(xmin', gridSize, 1) + repmat((xmax - xmin)', gridSize, 1) .* rand(gridSize, d);
Xgrid = [Xgrid; xx];
sx = size(Xgrid, 1);
maxes = zeros(1, nK);
[meanVector, varVector] = mean_var(Xgrid, xx, yy, KernelMatrixInv{1}, l(1,:), sigma, sigma0);

varVector(varVector<=sigma0 + eps) = sigma0 + eps;
sigVector = sqrt(varVector);

% Find the sample range [left, right].
probf = @(m0) prod(normcdf((m0 - meanVector)./sigVector));
left = max(meanVector);
leftprob = probf(left);
while leftprob > 0.1
    if left > 0.01
        left = left/2;
    else
        left = left*2 - 0.05;
    end
    leftprob = probf(left);
end

right = max(meanVector+5*sigVector);
rightprob = probf(right);
while (rightprob < 0.95)
    right = right + right - left;
    rightprob = probf(right);
end
mgrid = linspace(left, right, 100);


prob = prod(normcdf((repmat(mgrid,[sx,1]) - repmat(meanVector,[1,100]))./repmat(sigVector,[1, 100])),1);

if  isempty(find((prob>0.05) & (prob < 0.95), 1))
    maxes(1,:) = max(meanVector) + rand(1,nK)*0.1;
else
    if nK > 10
        % Randomly sample the function upper bounds from a Gumbel distribution
        % that approximates the CDF.
        med = find_between(0.5, probf, prob, mgrid, 0.1);
        
        q1 = find_between(0.25, probf, prob, mgrid, 0.1);
        
        q2 = find_between(0.75, probf, prob, mgrid, 0.1);
        beta=(q1-q2)/(log(log(4/3)) - log(log(4)));
        alpha = med+beta*log(log(2));
        assert(beta > 0);
        maxes(1,:) = - log( -log(rand(1, nK)) ) .* beta + alpha;
    else
        % Randomly sample the function upper bounds from the CDF.
        for k = 1:nK
            samp_rd = rand;
            if samp_rd > rightprob
                maxes(1,k) = right;
            elseif samp_rd < leftprob
                maxes(1,k) = left;
            else
                maxes(1,k) = find_between(samp_rd, probf, prob, mgrid, 0.1);
            end
        end
    end
    
end
maxmean  =  max(meanVector);
badidx = maxes(1,:) < maxmean + 5*sigma0;
maxes(1, badidx) = maxmean + 5*sigma0;

% Compute the acquisition function values on Xgrid.
gamma = (repmat(maxes(1,:),[sx 1]) - repmat(meanVector, [1, nK])) ./ repmat(sigVector, [1, nK]);

pdfgamma = normpdf(gamma);
cdfgamma = normcdf(gamma);
yvals =  sum(gamma.*pdfgamma./(2*cdfgamma) - log(cdfgamma),2);

[maxVal, maxIdx ] = max(yvals);
start = Xgrid(maxIdx,:);

% Optimize the acquisition function.
acfun = @(x) evaluateMES(x, maxes, xx, yy, KernelMatrixInv, l, sigma, sigma0);
neg_acfun = @(x) negative_wrapper(acfun, x);
[ optimum, fval] = fmincon(neg_acfun, start, [], [], [], [], xmin, xmax, [], ...
    optimset('MaxFunEvals', 100, 'TolX', eps, 'Display', 'off', 'GradObj', 'on'));
fval = -fval;
if fval < maxVal
    %disp('fmincon for MES-G failed to return a value better than the initialization.')
    [optimum, fval] = fminsearch(neg_acfun, start);
    fval = -fval;
    if fval < maxVal
        %disp('fminsearch for MES-G failed to return a value better than the initialization.')
        optimum = start;
    else
        % check if optimum is in range
        flag = 0;
        for i = 1:length(xmin)
            if optimum(i) < xmin(i) || optimum(i) > xmax(i)
                flag = 1;
                break;
            end
        end
        if flag
            %disp('fminsearch for MES-G failed to search within [xmin, xmax].')
            optimum = start;
        end
    end
end
