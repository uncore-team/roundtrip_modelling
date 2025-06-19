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
%    if (min(x) <= offset) % that model cannot assess these data
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


	% Based on D'Agostino p. 122

    logx = log(x - offset); % still ordered, now unshifted and normal
    w = (logx - mu) / sigma; % now normal(0,1)
    Z = normcdf(w,0,1); 

    [kindstat,stat] = StatisticGof(Z);
    if strcmp(kindstat,'A2')

        if modelnotfromdata
            % no correction of stat in this case (see table 4.2 D'Agostino)
    	    thresh = 2.492; % known parameters = previous model (n>=5)
                            % table 4.2 D'Agostino (confirmed by MonteCarlo in
                            % test_tabulategofthrs.m)
        else % model comes from data
            % correction: A2 for case 3 (both parameters were deduced from the same sample). 
            % the following is D'Agostino table 4.7, upper tail
            stat = stat * (1 + 0.75/n + 2.25/n^2);
            thresh = 0.752; 
        end

    elseif strcmp(kindstat,'W2')

        if modelnotfromdata
            % correction if parms are true (not from sample); table 4.2
            % D'Agostino
            stat = (stat - 0.4/n + 0.6/n^2) * (1 + 1/n);
            thresh = 0.461;  % we have confirmed this value with MonteCarlo (test_tabulategofthrs.m)
        else % model comes from data
            % the following is D'Agostino table 4.7, upper tail
            stat = stat * (1 + 0.5/n);

            % This is from D'Agostino, but does not work well for us since
            % D'Agostino reports for the Normal, not for the Lognormal.
            % thresh = 0.117;

            % This is the rsult of MonteCarlo of test_tabulate... :

            coeffs1 = [0.000000000000291  -0.000000000525001   0.000000372697872  -0.000130474298915   0.029244137137401   0.361152490347344]; % a 5th deg pol; get a norm of residuals = 0.59
            coeffs2 = [-0.000000000036426   0.000000928221916   0.018878529322011  -5.955872518883343]; % cubic fitting with normresid = 29.65
            numparts = 2;
            transweight = @(x,k,transitionpoint) 1 ./ (1 + exp(-k * (x - transitionpoint)));
            ks = 0.05 * ones(1,numparts - 1); % transition weights from one part to the next
            endspartsss = [600,10000];
            coeffsparts = {coeffs1,coeffs2};
            if n <= endspartsss(1)/2 
                parts = [1]; % polynomials to weld
            elseif n <= (endspartsss(1) + endspartsss(2))/2
                parts = [1,2];
                transitionpoint = endspartsss(1);
            else
                parts = [2];
            end
            
            if length(parts) == 1
                y(f) = polyval(coeffsparts{parts(1)},n);
            else
                poly1 = polyval(coeffsparts{parts(1)},n);
                poly2 = polyval(coeffsparts{parts(2)},n);     
                w = transweight(n,ks(parts(1)),transitionpoint);
                thresh = poly1 * (1 - w) + poly2 * w;
            end

        end

    else
        error('Unknown gof statistic');
    end

    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis of the data coming from an exponential
    end 

end
