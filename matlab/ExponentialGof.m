function [reject,stat,thresh] = ExponentialGof(x,alpha,beta)
% Based on D'Agostino p. 141: both parameters unknown. Same corrections as
% explained in ExponentialFit() apply here.

    if beta <= 0
        error('Invalid beta for exponential distr.');
    end

    n = length(x);
    %xsorted = sort(x);
    mu = 1/beta; % they use beta in the book when actually they want to use the mean 
    wsorted = (x - alpha) / mu;
    z = 1 - exp(-wsorted);
    Z = sort(z);

    % calculate statistic: A2 for case 3 (both parameters were deduced from
    % the same sample). This statistic measures the squared distance
    % between experimental and theoretical Zs, and, indirectly, between 
    % theoretical and experimental Xs (p. 100)
    sumatoria = sum(([1:n]*2-1).*log(Z)+(2*n+1-2*[1:n]).*log(1-Z)); % page 101, bottom formula
    A2 = -n - (1/n)* sumatoria;
    % do the following since parameters come from sample (D'Agostino table 4.14)
    A2 = A2*(1 + 5.4/n - 11/n^2); 

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
    
    % test the hypothesis 
    thresh = 1.321; % D'Agostino table 4.14
    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis of the data coming from an exponential
    end 

end