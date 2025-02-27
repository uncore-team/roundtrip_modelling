function cdf = LoglogisticCdf(a,b,c,x)
% Calculate the CDF of the loglogistic for the vector of values X

    cdf=1./(1+((x-a)/b).^(-1/c));
    cdf(x<=a)=0;

end