function [reject,stat,thresh] = LognormalGof(x,offset,mu,sigma,modelnotfromdata)

    LognormalCheckParms(offset,mu,sigma);

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

% % Taken from https://es.mathworks.com/matlabcentral/fileexchange/60147-normality-test-package
% % Paper of 2017 in tyrell project /ejecucion/docs/A Compilation of Some Popular Goodness of Fit Tests for Normal Distribution.pdf
% % They say they take it from D'Agostino p. 122 and table 4.9 in p.127
% 
%   Note that this only works for parameters coming from data, apparently.
%
%    if (min(x) < offset) % that model cannot assess these data
%        reject = 1;
%        stat = NaN;
%        thresh = NaN;
%        pvalue = 0;
%        return;
%    end
%
%     y = sort(log((x(:).') - offset)); % transform from lognormal to normal with expectation MU and std SIGMA
%     ui = normcdf(y,mu,sigma); % if taking mu,sigma from the data:  ui=normcdf(zscore(y),0,1); % zscore embed the mean and sigma estimated from the data
%     oneminusui = sort(1-ui);
%     i = 1:n;
%     lastt = (2*i-1).*(log(ui)+log(oneminusui));
%     asquare = -n-(1/n)*sum(lastt);
%     adj = 1+0.75/n+2.25/(n^2);
%     AD = asquare*adj;
% 
%     % this is from table 4.9 in D'Agostino and serves to calculate p-value, then compare it to the significance level
%     if AD<=0.2
%         pvalue=1-exp(-13.436+101.14*AD-223.73*AD^2);
%     elseif AD<=0.34
%         pvalue=1-exp(-8.318+42.796*AD-59.938*AD^2);
%     elseif AD<=0.6
%         pvalue=exp(0.9177-4.279*AD-1.38*AD^2);
%     elseif AD<=153.467
%         pvalue=exp(1.2937*AD-5.709*AD+0.0186*AD^2);
%     else
%         pvalue=0;
%     end
% 
%     stat = AD;
%     thresh = NaN;
%     if pvalue <= 0.05 % equivalently, the statistic is greater than the threshold
%         reject = 1; % the null hypothesis (data come from normal) should be rejected: the possible normality is caused by noise
%     else
%         reject = 0; % cannot reject the hypothesis of normality
%     end


	% Based on D'Agostino p. 122: parameters unknown.

    xord = sort(x - offset); % go to a non-shifted sample from a non-shifted lognormal
    if (xord(1) < 0) % that model cannot assess these data
                       % AD test cannot work with samples that produce 0 or 1 Z values
        reject = 1;
        stat = Inf;
        thresh = NaN;
        return;
    end
    logxord = log(xord); % still ordered, now normal
    m = mean(x);
    w = (logxord - mu) / sigma;
    Z = normcdf(w,0,1);
    A2 = ADstatistic(Z);
    if isnan(A2) % that model cannot assess these data
        reject = 1;
        stat = Inf;
        thresh = NaN;
        return;
    end
    % do the following if parameters come from sample (D'Agostino table 4.7)
    if ~modelnotfromdata
    	A2 = A2 * (1 + 0.75/n + 2.25/n^2);
    end
    stat = A2;
    
    % --- hypothesis test
    if modelnotfromdata
    	thresh = 2.492; % known parameters = previous model (n>=5)
                        % table 4.2 D'Agostino (confirmed by MonteCarlo in
                        % test_tabulategofthrs.m)
    else % unknown parameters, estimated from the very sample

        % the following is D'Agostino but does not work for use since he
        % assumes only 2 parameters unknown, but we have 3; thus we have
        % conducted new MonteCarlo experiments to deduce the threshold in
        % this case
        %
        % thresh = 0.752; % table 4.7, upper tail.

        % NOTE: We have tested with test_tabulategofthrs.m that, if we
        % assume the offset known and the rest of parms taken from the
        % sample, the threshold follows a straight line on the samplesize,
        % not a constant: 0.3123 * samplesize + 7.4905


        % MonteCarlo results in test_tabulategofthrs.m; fitting at the end
        % section of that script:

        % for the tabulated result using only samplesizes in [20,500]:
        %thresh = 1.177616088094618;

        % for the tabulated result using samplesizes in [20,1000]:
        thresh = 1.142624318695783;
        
    end
    if (stat > thresh) 
        reject=1; % reject
    else
        reject=0; % cannot reject
    end
    
end
