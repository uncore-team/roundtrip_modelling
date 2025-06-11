function cdf = LognormalCdf(offset,mu,sigma,x)
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offset,mu,sigma);

    cdf = logncdf(x-offset,mu,sigma);

end