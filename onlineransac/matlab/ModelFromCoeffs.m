function m = ModelFromCoeffs(v)
% Given a row vector V created with ModelToCoeffs, recover the model.
%
% V -> row vector with coeffs (see ModelToCoeffs).
% 
% M <- model (see ModelCreate).

    switch v(1)
        case 0 % LL3

            m = ModelCreate('LL3');
            m.coeffs.a = v(2);
            m.coeffs.b = v(3);
            m.coeffs.c = v(4);
            if (~isnan(m.coeffs.a)) && (~isnan(m.coeffs.b)) && (~isnan(m.coeffs.c)) && ...
               (~isinf(m.coeffs.a)) && (~isinf(m.coeffs.b)) && (~isinf(m.coeffs.c))
                m.defined = 1;
            end

        case 1 % EXP2

            m = ModelCreate('EXP2');
            m.coeffs.alpha = v(2);
            m.coeffs.beta = v(3);
            if (~isnan(m.coeffs.alpha)) && (~isnan(m.coeffs.beta)) && ...
               (~isinf(m.coeffs.alpha)) && (~isinf(m.coeffs.beta))
                m.defined = 1;
            end

        case 2 % LN3

            m = ModelCreate('LN3');
            m.coeffs.gamma = v(2);
            m.coeffs.mu = v(3);
            m.coeffs.sigma = v(4);
            if (~isnan(m.coeffs.gamma)) && (~isnan(m.coeffs.mu)) && (~isnan(m.coeffs.sigma)) && ...
               (~isinf(m.coeffs.gamma)) && (m.coeffs.sigma > 0)
                m.defined = 1;
            end

        case 3 % BERN

            m = ModelCreate('BERN');
            m.coeffs.ind0 = v(2);
            m.coeffs.ind1 = v(3);
            if (~isnan(m.coeffs.ind0)) && (~isnan(m.coeffs.ind1)) && ...
               (m.coeffs.ind1 >= m.coeffs.ind0)
                m.defined = 1;
            end

        otherwise
            error('Invalid model type');
    end

end