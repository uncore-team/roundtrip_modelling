function m = ModelCreate(ty)
% Create an undefined model with the given type.
%
% TY -> the type: 'LL3', 'EXP2', 'LN3', 'BERN'
%
% M <- the model. A struct:
%       .type <- copy of TY
%       .defined <- 1 if defined, 0 if not (0 after this function)
%       .coeffs <- a struct with the coefficients of the model.
%           --- If LL3:
%           .a,.b,.c
%           --- If EXP2:
%           .alpha,.beta
%           --- If LN3:
%           .gamma,.mu,.sigma
%           --- If BERN:
%           .ind0,.ind1

    m = struct('type',ty,...
               'defined',0);
    if strcmp(ty,'LL3')
        m.coeffs = struct('a',NaN,'b',NaN,'c',NaN);
    elseif strcmp(ty,'EXP2')
        m.coeffs = struct('alpha',NaN,'beta',NaN);
    elseif strcmp(ty,'LN3')
        m.coeffs = struct('gamma',NaN,'mu',NaN,'sigma',NaN);
    elseif strcmp(ty,'BERN')
        m.coeffs = struct('ind0',NaN,'ind1',NaN);
    else
        error('Unknown model');
    end

end