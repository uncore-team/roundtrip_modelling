function [reject,stat,thresh] = ModelGof(m,data)
% Check the goodness of fit of DATA w.r.t. the given model M.
%
% M -> model (see ModelFit()).
% DATA -> vector of data that have produced that model.
%
% REJECT <- 1 if we have to reject the null hypothesis that the model
%           explains the sample, 0 if we cannot say if it explain the sample
%			or not (and therefore we might assume the model is correct).
% STAT <- value of the statistic for that sample (the smaller, the better)
% THRESH <- value below which the statistic should have been in order not
%           to reject the null hypothesis with the 5% significance level, i.e.,
%			in order to assume the model fits well the data. In case of non-
%			rejection, STAT will be in [0,THRESH], thus it can be normalized to
%			[0,1] in order to get a "degree of non-rejection" (less rejection as
%			smaller is that degree); this works as an alternative to the p-value
%			that is safe to use as long as we do not assume any linearity or
%			other particular decreasing profile in STAT.

    if ~m.defined
        error('Undefined model cannot do gof');
    end

    if strcmp(m.type,'LL3')

        [reject,stat,thresh] = LoglogisticGoF(data,m.coeffs.a,m.coeffs.b,m.coeffs.c,0);

    elseif strcmp(m.type,'EXP2')

        [reject,stat,thresh] = ExponentialGof(data,m.coeffs.alpha,m.coeffs.beta);

    elseif strcmp(m.type,'LN3')

        [reject,stat,thresh,~] = LognormalGof(data,m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);

    elseif strcmp(m.type,'BERN')

        reject = 0;
        stat = NaN;
        thresh = NaN;

    else
        error('Invalid model type');
    end

end