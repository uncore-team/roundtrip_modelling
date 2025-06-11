function LognormalCheckParms(offs,mu,sigma)
% Produce an error if the parms are invalid.

    if (offs < 0) || (sigma <= 0)
        error('Invalid parameters for a lognormal');
    end

end