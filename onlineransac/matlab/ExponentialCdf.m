function cdf = ExponentialCdf(alpha,beta,xs)
% Calculate the CDF of the exponential for the vector of values X
% See the parameters in ExponentialFit.

	ExponentialCheckParms(alpha,beta);
    cdf = 1 - exp(-(xs - alpha)/beta);

end
