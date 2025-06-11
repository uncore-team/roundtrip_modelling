function v = LognormalToVariance(offset,mu,sigma)
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offset,mu,sigma);

    v = (exp(sigma^2) - 1) * exp(2 * mu + sigma^2);

end
