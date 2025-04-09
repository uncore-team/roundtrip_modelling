function m = ModelCreateRnd(type,mode)
% Create and define a model of type TYPE in mode MODE.
%
% MODE -> one of the following:
%           'typrnd' -> create the model with params drawn uniformly from
%                       the typical lower and upper bounds given by
%                       ModelTypParmBounds.

    if strcmp(mode,'typrnd')
        m = ModelCreate(type);
        ps = ModelTypParmBounds(type);
        nparms = size(ps,1);
        m.defined = 1; % just to force to gather a vector of coeffs
        coeffs = ModelToCoeffs(m);
        for f = 1:nparms
            coeffs(f + 1) = myrnd(ps(f,1),ps(f,2));
        end
        m = ModelFromCoeffs(coeffs);
    else
        error('Unknown model type');
    end

end

function num = myrnd(xmin, xmax)
    num = xmin + (xmax - xmin) * rand();
end 