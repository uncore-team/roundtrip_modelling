function v = ModelToVariance(m,scenario)
% Return the variance of that model.
%
% M -> model (see ModelCreate)
% SCENARIO -> complete scenario, needed for some models.
%
% V <- variance

    if ~m.defined
        error('Undefined model cannot have variance');
    end

    if strcmp(m.type,'LL3')

        v = LoglogisticToVariance(m.coeffs.a,m.coeffs.b,m.coeffs.c);

    elseif strcmp(m.type,'EXP2')

        v = ExponentialToVariance(m.coeffs.alpha,m.coeffs.beta);

    elseif strcmp(m.type,'LN3')

        v = LognormalToVariance(m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);

    elseif strcmp(m.type,'BERN')

        S = scenario(m.coeffs.ind0:m.coeffs.ind1);
        v = var(S);
        
    else
        error('Invalid model type');
    end  
end