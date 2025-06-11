function [reject,stat,thresh] = ExponentialGof(x,alpha,beta,modelnotfromdata)
% Based on D'Agostino p. 141: both parameters unknown. 
% See ExponentialFit()

    ExponentialCheckParms(alpha,beta);

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

	% ---- calculate the experimental EDF

    xsorted = sort(x);
    wsorted = (xsorted - alpha) / beta;
    Z = 1 - exp(-wsorted);

    % ---- calculate A2 statistic
    A2 = ADstatistic(Z);
    if isnan(A2) % that model cannot assess these data
        reject = 1;
        stat = Inf;
        thresh = NaN;
        return;
    end
    % correction: A2 for case 3 (both parameters were deduced from the same sample). 
    % do the following since parameters come from sample (D'Agostino table 4.14)
    if ~modelnotfromdata
        A2 = A2 * (1 + 5.4/n - 11/n^2); 
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
    	thresh = 2.492; % known parameters (n>=5)
                        % table 4.2 D'Agostino; confirmed with MonteCarlo
                        % experiments (test_tabulategofthrs.m)
    else % parameters come from the very data

        % the following is D'Agostino and it does work even when it 
        % assumes only one parameter unknown while we have 2; we have
        % conducted new MonteCarlo experiments to deduce the threshold in
        % this case in test_tabulategofthrs.m and got the same threshold
        % for sample sizes ranging from 20 to 10000.
        
        thresh = 1.321; 
        
    end    
    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis of the data coming from an exponential
    end 

end
