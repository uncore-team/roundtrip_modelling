function [alpha,beta,ok] = ExponentialFit(x)
% Estimates a shifted exponential that fits X in the MLE sense.
% This uses the shifted exponential defined in D'Agostino p. 133, from
% where we have: pdf(x;alpha,beta) = 1/beta * exp(-(x-alpha)/beta).
% expectation == beta + alpha; median == ln2 * beta + alpha; variance == beta^2
%
% Such an exponential distribution may generate samples with values equal
% to the offset.
%
% Wikipedia: pdf(x;alpha,lambda) = lambda * exp(-lambda*(x-alpha)); expectation == 1/lambda + alpha; median == ln2/lambda + alpha == ln2 * (expectation - alpha) + alpha
% Matlab: it has only the non-shifted distrib. where pdf(x,mu) = 1/mu * exp(-x/mu)
%
% From Wikipedia to Here: beta = 1/lambda
% From Here to Wikipedia: lambda = 1/beta
% From Here to Matlab: mu = beta
% From Matlab to Here: beta = mu (alpha = 0)
%
% This method always returns parameters, except when length(x) < certain
% minimum or some value in the sample is <= 0. In those cases alpha and 
% beta get NaN, and OK gets 0.
%
% X -> (possibly unordered) sample.
%
% ALPHA <- offset (location); it will be >= 0.
% BETA <- the non-shifted mean (shape); it will be > 0.
% OK <- 1 if there are enough values in the sample to do the fitting.

global TOLROUNDTRIPS

    ConstantsInit();

    ok = 0;
    alpha = NaN;
    beta = NaN;

    n = length(x);
    xord = sort(x);
    minx = xord(1);
    if minx <= 0
        warning('Tried to fit an exponential with min(x) <= 0');
        return;
    end

    % --- unbiased MLE estimation:
    mu = mean(x);
    beta = (mu - minx) * n / (n - 1); % estimate of the (non-shifted) mean
    alpha = minx - beta / n; % estimate of the offset
    if alpha < 0 % in that case we revert to the MLE, biased estimators
        % COMMENTED OUT because using the median breaks the MLE estimation
        % (it is no longer a MLE, i.e., does not maximize the likelihood 
        % any longer):
        %
        % % estimate beta through the median to be more robust under the very
        % % likely situation here of having large values in the right tail of
        % % the sample
        % if n/2 == floor(n/2) % xord has an even number of elements
        %     med = (xord(n/2) + xord(n/2 + 1)) / 2;
        % else % xord has an odd number of elements
        %     med = xord((n+1)/2);
        % end
        % beta = med / log(2);

        beta = mu - minx;
        alpha = minx - TOLROUNDTRIPS; % offset; slightly smaller than min(x) to avoid singularities in the AD test (we introduce a slight bias by doing this, but it is acceptable)    
        if alpha < 0
            alpha = 0; % should not happen since min(x) > 0
        end
    end

    % % ---- this is the method in D'Agostino p. 425 to reduce to a 0-offset exp.
    % % we do not use this since it is a biased estimation of both beta and
    % % alpha, as explained in docs/MLE estimation of EXP2.png
    % if n < 2
    %     waning('Tried to fit an exponential with too few values');
    %     return;
    % end
    % xordshifted = xord - minx; % then we get rid of the first value:
    % xordshifted = xordshifted(2:end); 
    % % the resulting sample is exponential with same shape but location 0;
    % % we use this new sample to estimate the shape. For that,
    % % use the median to reach the mean, since the median is a more robust
    % % estimator, i.e., we are using an estimator of the mean based on the 
    % % median
    % if n/2 == floor(n/2) % xordshifted has an odd number of elements
    %     med = xordshifted(n/2);
    % else % xordshifted has an even number of elements
    %     med = (xordshifted((n-1)/2) + xordshifted((n+1)/2)) / 2;
    % end
    % beta = med / log(2);
    % % then we go back to the original sample to estimate the offset:
    % % using the minimum of the sample is equivalent to the MLE estimation
    % % of the offset
    % alpha = minx - TOLROUNDTRIPS; % offset; slightly smaller than min(x) to avoid singularities in the AD test (we introduce a slight bias by doing this, but it is acceptable)    
    

    % % ---- this is my heuristic, and it is not formally correct:
    % 
    % alpha = xord - TOLROUNDTRIPS; % offset; slightly smaller than min(x) to avoid singularities in the AD test (we introduce a slight bias by doing this, but it is acceptable)    
    % xshifted = xord - alpha;
    % med = median(xshifted); % median(exponential) = (ln 2) / beta
    % % robust estimation because of the possible existance of too large values on the right tail that distort the mean:
    % beta = med/log(2);

    % % ---- This is an idea of selecting mean or median depending on the
    % % reliability of the sample, but i dimissed the idea:
    %
    % alpha = xord - TOLROUNDTRIPS; % offset; slightly smaller than min(x) to avoid singularities in the AD test (we introduce a slight bias by doing this, but it is acceptable)    
    % mea = mean(xshifted); % mean(exponential) = 1 / beta
    % % in theory, mean/median = 1/ln2 = 1.442...
    % 
    % if mea/med > 2 % we include some margin on the 1/ln(2) theoretical value
    %     % robust estimation because of probably too large values on the right:
    %     beta = med/log(2);
    % else
    %     % non-robust estimation:
    %     beta = mean(x - alpha); 
    % end


    % ---- Method according to D'Agostino, p. 141 
    % This method corrects the biases of both alpha and beta estimates, as
    % developed in docs/derivation of the unbiasing of EXP2 estimates.pdf

    % PROBLEM: This method produces negative alpha in some of
    % the estimations due to the fact that when the average of the sample
    % is > n * x1 then alpha is negative in this formulation; that can
    % occur due to large values in the sample, far to the right, that
    % distort the mean. 
    % It could be addressed by first estimating the
    % median and then doing mean = median / ln2, but that does not provide
    % any analytical guaranteee that that alpha < 0 in all cases (the
    % median might also be distorted).
    % We adddress the problem instead by falling back to the biased MLE
    % estimation in the cases that that occurs.

    ok = 1;

end