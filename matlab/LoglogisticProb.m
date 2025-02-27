function p = LoglogisticProb(a,b,c,x0,x1)
% Return the probability that the given LL3 lies within [x0,x1]

    if (x0>x1)
        error('Bad range for probability');
    end
    if (x0==x1)
        p=0;
        return;
    end
    if (x0==-Inf)
        x0=a;
    end
    if (x1<Inf)     
        cdf=LoglogisticCdf(a,b,c,[x0 x1]);
        p=cdf(2)-cdf(1);
    else
        cdf=LoglogisticCdf(a,b,c,x0);
        p=1-cdf(1);
    end
    
end
