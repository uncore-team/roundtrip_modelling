function cdf = LognormalCdf(offset,mu,sigma,x)
% Calculate the CDF of the lognormal for the vector of values X

    cdf = logncdf(x-offset,mu,sigma);

end