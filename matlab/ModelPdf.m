function ys = ModelPdf(m,scenario,xs)
% Calculate the pdf of the given model for the support in XS.
%
% M -> model (see ModelCreate()).
% SCENARIO -> complete scenario, needed for some models.
% XS -> vector of points in the support to calculate the pdf for.
%
% YS <- calculated pdf values for that support.

    if ~m.defined
        error('Undefined model cannot do pdf');
    end

    if strcmp(m.type,'LL3')

        ys = LoglogisticPdf(xs,m.coeffs.a,m.coeffs.b,m.coeffs.c);

    elseif strcmp(m.type,'EXP2')

        ys = ExponentialPdf(xs,m.coeffs.alpha,m.coeffs.beta);

    elseif strcmp(m.type,'LN3')

        ys = LognormalPdf(xs,m.coeffs.gamma,m.coeffs.mu,m.coeffs.sigma);

    elseif strcmp(m.type,'BERN')

        S = scenario(m.coeffs.ind0:m.coeffs.ind1);
        [ys,~] = histcounts(S,[xs - min(diff(xs))/2, xs(end) + min(diff(xs))/2],'Normalization','pdf');        
        
    else
        error('Invalid model type');
    end  

end