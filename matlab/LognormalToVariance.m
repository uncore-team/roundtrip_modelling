function v = LognormalToVariance(gamma,mu,sigma)
% Return the variance of that lognormal.

    v = (exp(sigma^2) - 1) * exp(2 * mu + sigma^2);

end
