function ds = ModelRnd(mo,n,m)
% Generate a matrix of NxM random numbers drawn from MO.
% The resulting sample may have values equal to the offset (if the model 
% has offset) due to 2 reasons: they are naturally equal to it, if the 
% model allows for that; they become equal to it when
% adding a large offset to a small value, due to numerical imprecissions.
% The latter can happens also in a real scenario when we are measuring delays:
% some of them may come from a continuous support distrib., i.e., infinite
% precision, but become truncated in their precision.

    if ~mo.defined
        error('Cannot draw numbers from undefined model');
    end

    if strcmp(mo.type,'LL3')
        ds = LoglogisticRnd(mo.coeffs.a,mo.coeffs.b,mo.coeffs.c,n,m);
    elseif strcmp(mo.type,'LN3')
        ds = LognormalRnd(mo.coeffs.gamma,mo.coeffs.mu,mo.coeffs.sigma,n,m);
    elseif strcmp(mo.type,'EXP2')
        ds = ExponentialRnd(mo.coeffs.alpha,mo.coeffs.beta,n,m);
    else
        error('Unknown model type');
    end

end