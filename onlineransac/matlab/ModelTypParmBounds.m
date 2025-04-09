function bs = ModelTypParmBounds(t)
% Given in T a model type (see ModelCreate), return a matrix with as many
% rows as parameters with 2 columns: lower and upper bounds of each
% parameter in typical situations.

    if strcmp(t,'LL3')
        bs = [1e-4 76e3; 1e-4 32e3; 0.05 0.45];
    elseif strcmp(t,'LN3')
        bs = [1e-4 76e3; -10 10; 1e-3 8];
    elseif strcmp(t,'EXP2')
        bs = [1e-3 26e3; 1e-5 5];
    else
        error('Unknown model type');
    end

end