function r = LognormalRnd(offs, mu, sigma, m, n)
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offs,mu,sigma);
    
    r = lognrnd(mu,sigma,m,n) + offs;

end
