function [E,V]=LognormalToEV(offs,mu,sigma)
% Return the expectation and variance for the 3-logn, which will include the
% contribution of the offset.
% See lognormal parms in LognormalFit.m

    LognormalCheckParms(offs,mu,sigma);

    E = LognormalToExpectation(offs,mu,sigma);
    V = LognormalToVariance(offs,mu,sigma);

end
