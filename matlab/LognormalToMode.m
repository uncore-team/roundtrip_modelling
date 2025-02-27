function m = LognormalToMode(offset,mu,sigma)
% return the mode of that 3-lognormal

    if (sigma <= 0)
        error('Invalid lognormal');
    end
    m = exp(mu - sigma^2) + offset;

end