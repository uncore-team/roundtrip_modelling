function ex = LognormalToExpectation(gamma,mu,sigma)
% Return the expectation of that lognormal

    ex = gamma + exp(mu + (sigma^2)/2);

end