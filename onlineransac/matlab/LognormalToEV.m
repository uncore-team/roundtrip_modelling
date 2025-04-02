function [E,V]=LognormalToEV(gamma,mu,sigma)
% Return the expectation and variance for the 3-logn, which will include the
% contribution of the offset.

    E = gamma + exp(mu + (sigma^2)/2);
    V = (exp(sigma^2) - 1) * exp(2 * mu + sigma^2);

end
