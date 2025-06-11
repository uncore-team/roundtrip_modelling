function ex = LognormalToExpectation(offs,mu,sigma)
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offs,mu,sigma);

    ex = offs + exp(mu + (sigma^2)/2);

end