% This function is from the code for the paper
% Hernández-Lobato J. M., Hoffman M. W. and Ghahramani Z.
% Predictive Entropy Search for Efficient Global Optimization of Black-box
% Functions, In NIPS, 2014.
% See also: posteriorMean.m, gradientPosteriorMean.m at
% https://bitbucket.org/jmh233/codepesnips2014
function [f, g] = posteriorMean(x, Xsamples, Ysamples, KernelMatrixInv, l, sigma)
% Compute the posterior mean function
f = 0;
for i = 1 : size(KernelMatrixInv, 2)
    kstar = computeKnm(x, Xsamples, l(i,:)', sigma(i));
    f = f + kstar * KernelMatrixInv{ i } * Ysamples;
end
f = f / size(KernelMatrixInv, 2);
if nargout >= 2
    % Compute the gradient of the posterior mean function
    g = 0;
    for i = 1 :  size(KernelMatrixInv, 2)
        kstar = computeKnm(x, Xsamples, l(i,:)', sigma(i));
        dkstar = -repmat(l(i,:), size(Xsamples, 1), 1) .* ...
            (repmat(x, size(Xsamples, 1), 1) - Xsamples) .* repmat(kstar', 1, size(Xsamples, 2));
        g = g + dkstar' * KernelMatrixInv{ i } * Ysamples;
    end
    g = g / size(KernelMatrixInv, 2);
    
end