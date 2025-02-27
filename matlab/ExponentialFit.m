function [alpha,beta] = ExponentialFit(x)
% According to D'Agostino, p. 141 but corrected for the final beta to be
% 1/mean (not correct in the book, where they say they estimate beta when 
% actually they are estimating the mean).
% Shifted exponential pdf in D'Agostino p. 133
% 1/beta == mean of the distribution
% alpha == location of the distribution

    minlen = 10;

    n = length(x);
	if (n < minlen)
		error('Cannot fit anything with less than %d values',minlen);
	end
    
    mi = min(x);
    mu = n * (mean(x) - mi) / (n - 1); % estimate of the (non-shifted) mean
    alpha = mi - mu / n;

    beta = 1/mu; % beta is the reciprocal of the mean (in the book they use beta as the mean)

end