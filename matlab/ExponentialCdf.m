function cdf = ExponentialCdf(alpha,beta,xs)
% Calculate the CDF of the exponential for the vector of values X

    cdf = 1 - exp(-beta * (xs - alpha));

end