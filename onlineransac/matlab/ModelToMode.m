function mo = ModelToMode(m,scenario)
% Return the mode of the model.
% SCENARIO -> complete scenario, needed for some models.
%
% M -> model (see ModelCreate)
%
% MO <- mode

    if ~m.defined
        error('Undefined model cannot have mode');
    end

    if strcmp(m.type,'LL3')

        mo = LoglogisticToMode(m.coeffs.a,m.coeffs.b,m.coeffs.c);

    elseif strcmp(m.type,'EXP2')

        mo = ExponentialToMode(m.coeffs.alpha,m.coeffs.beta);

    elseif strcmp(m.type,'LN3')

        mo = LognormalToMode(m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);

    elseif strcmp(m.type,'BERN')

        S = scenario(m.coeffs.ind0:m.coeffs.ind1);
        mo = mode(S);
        
    else
        error('Invalid model type');
    end  

end