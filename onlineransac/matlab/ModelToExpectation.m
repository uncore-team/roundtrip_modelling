function ex = ModelToExpectation(m,scenario)
% Return the expectation of the model
%
% M -> model (see ModelCreate)
% SCENARIO -> complete scenario, needed for some models.
%
% EX <- expectation

    if ~m.defined
        error('Undefined model cannot have expectation');
    end

    if strcmp(m.type,'LL3')

        ex = LoglogisticToExpectation(m.coeffs.a,m.coeffs.b,m.coeffs.c);

    elseif strcmp(m.type,'EXP2')

        ex = ExponentialToExpectation(m.coeffs.alpha,m.coeffs.beta);

    elseif strcmp(m.type,'LN3')

        ex = LognormalToExpectation(m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);

    elseif strcmp(m.type,'BERN')

        S = scenario(m.coeffs.ind0:m.coeffs.ind1);
        ex = mean(S);
        
    else
        error('Invalid model type');
    end  

end