function [alpha,beta] = ExponentialFit(x)
% Using a simple MLE estimation that forces alpha > 0

    alpha = min(x) - eps; % offset; slightly smaller than min(x) to avoid singularities in the AD test (we introduce a slight bias by doing this, but it is acceptable)
    beta = 1.0 / mean(x - alpha); % beta is the reciprocal of the mean (in the book they use beta as the mean)

% % According to D'Agostino, p. 141 but corrected for the final beta to be
% % 1/mean (not correct in the book, where they say they estimate beta when 
% % actually they are estimating the mean).
% % Shifted exponential pdf in D'Agostino p. 133
% % 1/beta == mean of the distribution (no offset)
% % alpha == location of the distribution
% 
% % PROBLEM: This method produces negative alpha in a high proportion of
% the estimations
% 
%     minlen = 10;
% 
%     n = length(x);
% 	if (n < minlen)
% 		error('Cannot fit anything with less than %d values',minlen);
% 	end
% 
%     mi = min(x);
%     mu = n * (mean(x) - mi) / (n - 1); % estimate of the (non-shifted) mean
%     alpha = mi - mu / n; % offset
% 
%     beta = 1/mu; % beta is the reciprocal of the mean (in the book they use beta as the mean)

end