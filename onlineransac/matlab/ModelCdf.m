function ys = ModelCdf(m,scenario,xs)
% Calculate the cdf of the given model for the support in XS.
%
% M -> model (see ModelCreate()).
% SCENARIO -> complete scenario, needed for some models.
% XS -> vector of points in the support to calculate the cdf for.
%
% YS <- calculated cdf values for that support.

    if ~m.defined
        error('Undefined model cannot do pdf');
    end

    if strcmp(m.type,'LL3')

        ys = LoglogisticCdf(m.coeffs.a,m.coeffs.b,m.coeffs.c,xs);

    elseif strcmp(m.type,'EXP2')

        ys = ExponentialCdf(m.coeffs.alpha,m.coeffs.beta,xs);

    elseif strcmp(m.type,'LN3')

        ys = LognormalCdf(m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma,xs);

    elseif strcmp(m.type,'BERN')

        S = scenario(m.coeffs.ind0:m.coeffs.ind1);
        [pdf,~] = histcounts(S,[xs - min(diff(xs))/2, xs(end) + min(diff(xs))/2],'Normalization','pdf');        
        n = length(xs);
        ys = zeros(1,n);
        for f = 2:n
            ys(f) = ys(f-1) + (xs(f) - xs(f-1)) * (pdf(f-1) + (pdf(f) - pdf(f-1))/2); % trapezoidal approx
        end
        
    else
        error('Invalid model type');
    end  

end