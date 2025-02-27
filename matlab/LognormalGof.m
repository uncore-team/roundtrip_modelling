function [reject,stat,thresh,pvalue] = LognormalGof(x,offset,mu,sigma)

    n = length(x);
    if (sigma <= 0) 
        error('Invalid lognormal distr.');
    end
    if (min(x) < offset) % that model cannot assess these data
        reject = 1;
        stat = NaN;
        thresh = NaN;
        pvalue = 0;
        return;
    end

% Taken from https://es.mathworks.com/matlabcentral/fileexchange/60147-normality-test-package

    y = sort(log((x(:).') - offset)); % transform from lognormal to normal with expectation MU and std SIGMA
    ui = normcdf(y,mu,sigma); % if taking mu,sigma from the data:  ui=normcdf(zscore(y),0,1); % zscore embed de mean and sigma estimated from the data
    oneminusui = sort(1-ui);
    i = 1:n;
    lastt = (2*i-1).*(log(ui)+log(oneminusui));
    asquare = -n-(1/n)*sum(lastt);
    adj = 1+0.75/n+2.25/(n^2);
    AD = asquare*adj;
     
    if AD<=0.2
        pvalue=1-exp(-13.436+101.14*AD-223.73*AD^2);
    elseif AD<=0.34
        pvalue=1-exp(-8.318+42.796*AD-59.938*AD^2);
    elseif AD<=0.6
        pvalue=exp(0.9177-4.279*AD-1.38*AD^2);
    elseif AD<=153.467
        pvalue=exp(1.2937*AD-5.709*AD+0.0186*AD^2);
    else
        pvalue=0;
    end

    stat = AD;
    thresh = NaN;
    if pvalue <= 0.05 % equivalently, the statistic is greater than the threshold
        reject = 1; % the null hypothesis (data come from normal) should be rejected: the possible normality is caused by noise
    else
        reject = 0; % cannot reject the hypothesis of normality
    end

% Based on D'Agostino p. 122: parameters unknown.
%     xord = sort(x);     
%     logxord = log(xord); % still ordered, now normal
%     m = mean(x);
%     w = (logxord - mu) / sigma;
%     Z = normcdf(w,0,1);
% 
%     % calculate statistic: A2 for case 3 (both parameters were deduced from
%     % the same sample). This statistic measures the squared distance
%     % between experimental and theoretical Zs, and, indirectly, between 
%     % theoretical and experimental Xs (p. 100)
%     sumatoria = sum(([1:n]*2-1).*log(Z)+(2*n+1-2*[1:n]).*log(1-Z)); % page 101, bottom formula
%     A2 = -n - (1/n)* sumatoria;
%     % do the following since parameters come from sample (D'Agostino table 4.7)
%     A2 = A2*(1 + 0.75/n + 2.25/n^2); 
%     stat = A2;
% 
%     % test the hypothesis 
%     thresh = 2.492; % table 4.2 D'Agostino
%     if (stat > thresh) 
%         reject=1; % reject
%     else
%         reject=0; % cannot reject
%     end 

end
