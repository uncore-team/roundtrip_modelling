function mo = ModelChangeParms(m,parms)
% Given a model M, define it with the parms in PARMS.
%
% M -> a model as in ModelCreate.
% PARMS -> a vector with as many values as parameters.
%
% MO <- same as M but defined with those parameters.

    m2 = m;
    m2.defined = 1; % just to force dumping of parms
    coeffs = ModelToCoeffs(m2);
    coeffs(2:2 + length(parms) - 1) = parms;
    mo = ModelFromCoeffs(coeffs);

end
