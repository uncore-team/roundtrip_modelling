function [reject,stat,thresh] = ExponentialGof(x,alpha,beta,flagprevmodel)
% Based on D'Agostino p. 141: both parameters unknown. Same corrections as
% explained in ExponentialFit() apply here.

    % consider by default that alpha and beta have been estimated from
    % samples.
    if nargin == 3
        flagprevmodel = 0;
    end

    if alpha < 0
        error('Invalid alpha for exponential distr.: alpha < 0');
    end
    if beta <= 0
        error('Invalid beta for exponential distr.: beta <= 0');
    end

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

	% ---- calculate the experimental EDF

    xsorted = sort(x);
    if (xsorted(1) < alpha) % that model cannot assess these data
        reject = 1;
        stat = Inf;
        thresh = NaN;
        return;
    end
    mu = 1 / beta; % they use beta in the book when actually they want to use the mean 
    wsorted = (xsorted - alpha) / mu;
    z = 1 - exp(-wsorted);
    Z = sort(z);

    % --- calculate the A2 statistic: This statistic measures the squared distance
    % between experimental and theoretical Zs, and, indirectly, between 
    % theoretical and experimental Xs (p. 100)
    ind = 1:n;
    sumatoria = sum((2*ind-1).*log(Z)+(2*n+1-2*ind).*log(1-Z)); % page 101, bottom formula
    A2 = -n - (1/n)*sumatoria;
    % correction: A2 for case 3 (both parameters were deduced from the same sample). 
    % do the following since parameters come from sample (D'Agostino table 4.14)
    if ~flagprevmodel
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
    if flagprevmodel
    	thresh = 2.492; % known parameters = previous model (n>=5)
                        % table 4.2 D'Agostino
    else
        if n > 100
    	    thresh = 1.321; % long samples
        else
    	    ns = [5,10,15,20,25,50,100]; % lengths listed in D'Agostino table 4.15
    	    ts = [0.725,0.920,1.009,1.062,1.097,1.197,1.250];
    	    thresh = spline(ns,ts,n); % best interpolation since points are non-linear
    								    % do not include the point at > 100 since it is
    								    % too far and distorts the spline strongly
        end
    end    
    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis of the data coming from an exponential
    end 

end
