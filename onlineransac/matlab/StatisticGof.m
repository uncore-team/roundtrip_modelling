function [kind,stat] = StatisticGof(Z)
% Calculate an EDF statistic appropriate for the given Z (D'Agostino p.
% 97-101). The statistic measures the distance between experimental and 
% theoretical Zs, and, indirectly, between theoretical and experimental Xs.
%
% The statistic follows certain right-tailed distribution. We can set in
% that distribution a threshold value (in its support)
% corresponding to a given significance level. 
% Then, if the value calculated for the statistic falls 
% above that threshold, the goodness of fit hypothesis should be rejected
% (this is easier as the significance level grows).
% The p-value is the probability of the statistic distribution to
% produce a value of the statistic equal to or greater than the
% calculated one. The p-value will shrink as more strongly rejected
% is the null hypothesis. 
%
% The Anderson–Darling test tends to be more powerful than the Cramér–von 
% Mises test for most alternatives, especially when discrepancies occur in 
% the tails. Anderson-Darling gives more weight to the tails of the 
% distribution, while Cramér-von Mises treats all parts of the distribution 
% more equally. Therefore, for our distributions, it seems that AD has more
% power (higher prob. of rejecting the fit if the sample does not come from
% the required distribution). For instance, in https://apps.dtic.mil/sti/tr/pdf/ADA238630.pdf
% they test the power against uniform-produced samples and get more power
% for the AD.
%
% Z -> sample transformed to be theoretically uniformly distributed.
%
% KIND <- 'A2' for Anderson-Darling A2 statistic; used if Z does not
%              contain values too close to 0 or 1.
%         'W2' for the Cramer-von Mises statistic; used in the rest of
%              cases.
% STAT <- value of the statistic.

    n = length(Z);
    if n < 2
        error('Less than 2 Z values for calculating gof statistic');
    end
    Zsorted = sort(Z);
    is = 1:n;
    
    if (Zsorted(1) < eps) || (1 - Zsorted(end) < eps) 

        % W2 statistic (valid for samples that have 0 or 1 in Z)
        sumatoria = sum((Zsorted - (2*is - 1)/(2 * n)).^2); % p. 101, formula 4.2
        stat = sumatoria + 1/(12 * n);
        kind = 'W2';

    else 

        % A2 statistic (cannot work with samples that have 0 or 1 Z values)
        sumatoria = sum((is*2-1).*log(Zsorted)+(2*n+1-2*is).*log(1-Zsorted)); % page 101, bottom formula
        stat = -n - (1/n)* sumatoria;
        kind = 'A2';

    end

end
