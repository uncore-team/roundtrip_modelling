function [nps,nams] = ModelParmsDef(ty)
% Given a type of model, return the number of parameters (NPS) and their
% names in a cell (NAMS).

    if strcmp(ty,'LL3')
        nps = 3;
        nams = {'a','b','c'};
    elseif strcmp(ty,'LN3')
        nps = 3;
        nams = {'gamma','mu','sigma'};
    elseif strcmp(ty,'EXP2')
        nps = 2;
        nams = {'alpha','beta'};
    else
        error('Unknown model name');
    end

end