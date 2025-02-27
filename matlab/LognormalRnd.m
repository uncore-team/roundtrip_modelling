function r = LognormalRnd(offs, mu, sigma, m, n)
% generate m*n data from lognormal (offs,mu,sigma)

	r = lognrnd(mu,sigma,m,n) + offs;

end
