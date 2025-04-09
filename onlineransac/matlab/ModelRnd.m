function ds = ModelRnd(mo,n,m)
% Generate a matrix of NxM random numbers drawn from MO.

    if ~mo.defined
        error('Cannot draw numbers from undefined model');
    end

    if strcmp(mo.type,'LL3')
        ds = LoglogisticRnd(mo.coeffs.a,mo.coeffs.b,mo.coeffs.c,n,m);
    elseif strcmp(mo.type,'LN3')
        ds = LognormalRnd(mo.coeffs.gamma,mo.coeffs.mu,mo.coeff.sigma,n,m);
    elseif strcmp(mo.type,'EXP2')
        ds = ExponentialRnd(mo.coeffs.alpha,mo.coeffs.beta,n,m);
    else
        error('Unknown model type');
    end

end