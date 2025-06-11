function m = LognormalToMode(offset,mu,sigma)
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offset,mu,sigma);

    m = exp(mu - sigma^2) + offset;

end