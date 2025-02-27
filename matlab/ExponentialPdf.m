function ys = ExponentialPdf(xs,alpha,beta)
% Shifted exponential pdf in D'Agostino p. 133
% 1/beta == mean of the distribution
% alpha == location of the distribution
% pdf(x) == beta * exp(-beta * (x - alpha))
    
    if beta <= 0
        error('Invalid beta for exponential distr.');
    end
    
    ys = beta * exp(-beta * (xs - alpha));

end