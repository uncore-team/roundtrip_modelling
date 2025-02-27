function x = ExponentialRnd(alpha,beta,n,m)
% Shifted exponential pdf in D'Agostino p. 133
% 1/beta == mean of the distribution
% alpha == location of the distribution

    if beta <= 0
        error('Invalid beta for exponential distr.');
    end
    x = exprnd(1/beta,n,m) + alpha;

end