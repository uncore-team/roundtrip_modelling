function [reject,stat,thresh]=LoglogisticGoF(x,a,b,c,modelnotfromdata)
% Anderson-Darling test the goodness of fit of the 3-loglogistic (A,B,C) for 
% the sample XS, from which the very parameters have been deduced, using a 0.05
% significance level, according to "Goodnes-of-fit techniques", D'Agostino and 
% Stephens (1986)
%
% According to D'Agostino, gof tests do not seek to reject the null hypothesis
% (and thus to support the contrary hypothesis) because for that we should
% provide a null hypothesis of the style "all the distributions that are not LL"
% which is impossible. Therefore, gof tests are focused on non-rejections, given
% as null hypothesis the particular distribution we wish to test. If the test
% rejects it, we know that the distribution does not hold (under the significance
% level), and if the tests does not reject it, we cannot say that it does not
% hold, which is taken as enough support for the distribution to hold (but it could
% be that a different distribution holds strongerly).
%
% For the rest of the matters, a gof test is just like any other hypothesis 
% test: we define a certain statistic associated to the null hypothesis 
% (a statistic is a function of the observation data) for which we know a (usually
% right-tailed) distribution; if the probability, according to that distribution, of
% the statistics having the observed value or greater is small enough (e.g., 5%)
% then we reject the null hypothesis. Thus, the statistic becomes an index of
% how much the observations agree with the null hypothesis: the smaller the
% statistic (to the left of the support of the statistic distribution), the 
% greater the agreement. The p-value is the probability of observing
% that value of the statistic (i.e., those observation data) or greater, given 
% that the null hypothesis is true, or in other words, the area of the statistic
% distribution from the statistic value to the right); therefore, the smaller 
% the p-value, the smaller the support of the null hypothesis -> greater 
% statistics produce smaller p-values. 
% The p-value is calculated upon the distribution of the statistic.
%
% In a gof test, we are interested in large p-values (small statistics) in order
% not to reject the null hypothesis of the distribution for the data. In this gof
% in particular, the 5% p-value (significance level) corresponds to a statistic 
% of 0.660, but there is no actual calculation of the p-value: we just not-reject 
% the hypothesis if the statistic is less than that value, and reject it if it 
% is greater. 

% According to D'Agostino, p. 91, when there is a specific EDF-based gof test for
% a given distribution (which is the case for the LL), it is better to use it 
% than using a multi-purpose chi-squared type test, because we get more power.
% Actually, Anderson-Darling is more powerful than other EDF tests (D'Agostino
% p. 110) because it tests better the right tail of heavy-tailed distributions.
%
% REJECT <- 1 if we have to reject the null hypothesis that the model
%           explains the sample, 0 if we cannot say if it explain the sample
%			or not (and therefore we might assume the model is correct).
% STAT <- value of the statistic for that sample (the smaller, the better)
% THRESH <- value below which the statistic should have been in order not
%           to reject the null hypothesis with the 5% significance level, i.e.,
%			in order to assume the LL3 model fits well the data. In case of non-
%			rejection, STAT will be in [0,THRESH], thus it can be normalized to
%			[0,1] in order to get a "degree of non-rejection" (less rejection as
%			smaller is that degree); this works as an alternative to the p-value
%			that is safe to use as long as we do not assume any linearity or
%			other particular decreasing profile in STAT.
% modelnotfromdata -> 1 if the parameters do not come from sample; 0 if they have
%                  been calculated from the same sample.

    if (a < 0) 
        error('Invalid loglogistic distr.');
    end

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

    % ---- transform sample to LL2 (0 offset) and model (a,b,c) into Matlab model (mu,sigma)
    xLL = sort(x); % ordered needed later on
    if (xLL(1) <= a)
    	reject = 1; % cannot accept a distribution if some value falls below or at its minimum
    	stat = Inf;

    	thresh = NaN;
    	return;
    end
    xL = log(xLL - a); % since XLL is ordered and log does not change the order, xL is ordered
    mu = log(b); % converting from a,b,c to D'Agostino (alpha in the book), which uses the same as Matlab
    s = c; % converting from a,b,c to D'Agostino (beta in the book), which uses the same as Matlab

    % ---- calculate a new random variable Z (p. 160) that has uniform distribution in [0,1]
    % if the theoretical model (mu,sigma) is true (p. 101)
    % NOTE: Here, xL must be ordered 
    Z = 1./(1 + exp(-(xL-mu)./s)); % cdf formula for logistic

    % ---- calculate statistic: A2 
    A2 = statisticA2(Z);
    % if isnan(A2) % that model cannot assess these data
    %     reject = 1;
    %     stat = Inf;
    %     thresh = NaN;
    %     return;
    % end
    if ~modelnotfromdata
        % do the following only if parameters come from sample:
        A2 = A2*(1+.25/n);  % correction needed because both parameters are deduced from the same sample, table 4.22 case 3
    end
    stat = A2; % this statistic follows certain right-tailed distribution. We can set in
               % that distribution a threshold value (in its support)
               % corresponding to a given significance level. 
               % Then, if the value calculated for the statistic falls 
               % above that threshold, the hypothesis should be rejected
               % (this is easier as the significance level grows).
               % The p-value is the probability of the statistic distribution to
               % produce a value of the statistic equal to or greater than the
               % calculated one. The p-value will shrink as more strongly rejected
               % is the null hypothesis. We do not calculate it here
               % because the distribution of the statistic is not found in
               % the book.
    
    % ---- test the hypothesis 
    if modelnotfromdata
        thresh = 2.492; % for the case that parameters do not come from sample; 0.05, p. 105, table 4.2, right tail
                        % Confirmed by MonteCarlo (test_tabulategofthrs.m)
    else

        % the following is D'Agostino but does not work for use since he
        % assumes only 2 parameters unknown, but we have 3; thus we have
        % conducted new MonteCarlo experiments to deduce the threshold in
        % this case
        %
        % thresh = 0.660; % threshold for the case of parameters coming from sample; 0.05 significance level; page 157, table 4.22, case 3

        % MonteCarlo results in test_tabulategofthrs.m; fitting at the end
        % section of that script:

        if n < 90 % firs part: 5th poly

            parms = [-0.000000000421424;0.000000118010409;-0.000004990075374;-0.001491034210112;0.192700180267995;-0.707864114472733];
            thresh = parms(1) * n^5 + parms(2) * n^4 + parms(3) * n^3 + parms(4) * n^2 + parms(5) * n + parms(6);

        elseif n < 120 % link between first and second part: linear interpolation

            x1 = 90;
            y1 = 6.174206416983618; % value of the first part at 90
            x2 = 120;
            y2 = 6.237947960378145; % value of the second part at 120
            thresh = (y2 - y1)/(x2 - x1) * (n - x1) + y1;

        elseif n < 380 % second part: horizontal

            thresh = 6.237947960378145;

        else % third part: decaying linear

            m = -0.002007806860578; % slope of a line that starts at x = 0
            thresh = m * (n - 380) + 6.237947960378145; % line that starts at n = 380 (second part)

        end

    end 
    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis that the data come from a LL3
    end 

end
