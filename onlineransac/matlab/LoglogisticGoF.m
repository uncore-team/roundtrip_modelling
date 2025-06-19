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

    LoglogisticCheckParms(a,b,c);

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

    % ---- transform sample to LL2 (0 offset) and model (a,b,c) into Matlab model (mu,sigma)
    xL = log(x - a); % de-shift and de-exp
    mu = log(b); % converting from a,b,c to D'Agostino (alpha in the book), which uses the same as Matlab
    s = c; % converting from a,b,c to D'Agostino (beta in the book), which uses the same as Matlab

    % ---- calculate a new random variable Z (p. 160) that has uniform distribution in [0,1]
    % if the theoretical model (mu,sigma) is true (p. 101)
    Z = 1./(1 + exp(-(xL-mu)./s)); % cdf formula for logistic

    % ---- calculate statistic: A2 
    [kindstat,stat]  = StatisticGof(Z);
    if strcmp(kindstat,'A2')

        if modelnotfromdata
            % no correction of stat in this case (see table 4.2 D'Agostino)
            thresh = 2.492; % for the case that parameters do not come from sample; 0.05, p. 105, table 4.2, right tail
                            % Confirmed by MonteCarlo (test_tabulategofthrs.m)
        else % model comes from data
            % do the following only if parameters come from sample
            % correction needed because both parameters are deduced from the same sample, table 4.22 case 3
            stat = stat*(1+.25/n);  
            % the following is D'Agostino but does not work well; no
            % montecarlo has been done for this since we have turned to W2,
            % thus this will not work!
            thresh = 0.660; % threshold for the case of parameters coming from sample; 0.05 significance level; page 157, table 4.22, case 3

        end

    elseif strcmp(kindstat,'W2')

        if modelnotfromdata
            % correction if parms are true (not from sample); table 4.2
            % D'Agostino
            stat = (stat - 0.4/n + 0.6/n^2) * (1 + 1/n);
            thresh = 0.461; % we have confirmed this value with MonteCarlo (test_tabulategofthrs.m)
        else % model comes from data
            % the following is D'Agostino table 4.22
            stat = (n * stat - 0.08) / (n - 1);
            thresh = 0.098;
        end

    else
        error('Unknown gof statistic');
    end

    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis that the data come from a LL3
    end 

end
