function m = ModelFit(data,ind0,ind1,ty)
% Define a model for the data through fitting.
%
% DATA -> a vector of data.
% IND0, IND1 -> first and last indexes of the data within the complete 
%               scenario.
% TY -> type of model (see ModelCreate()).
%
% M <- model (see ModelCreate()). If fitting didn't work, the model is
%      undefined.

    if ind1 - ind0 + 1 ~= length(data)
        error('Invalid indexes of sample');
    end

    m = ModelCreate(ty);
    if strcmp(ty,'LL3')

        [a, b, c, exitflag] = LoglogisticFit(data);
        if (exitflag >= 0) 
            m.defined = 1;
            m.coeffs.a = a;
            m.coeffs.b = b;
            m.coeffs.c = c;
        end

    elseif strcmp(ty,'EXP2')

        [alpha,beta] = ExponentialFit(data);
        m.defined = 1;
        m.coeffs.alpha = alpha;
        m.coeffs.beta = beta;

    elseif strcmp(ty,'LN3')

        [ok, offs, mu, sigma] = LognormalFit(data);
        if ok
            m.defined = 1;
            m.coeffs.gamma = offs;
            m.coeffs.mu = mu;
            m.coeffs.sigma = sigma;
        end

    elseif strcmp(ty,'BERN')

        m.defined = 1;
        m.coeffs.ind0 = ind0;
        m.coeffs.ind1 = ind1; % Bernoulli model equals the very sample
        
    else
        error('Invalid model type');
    end

end