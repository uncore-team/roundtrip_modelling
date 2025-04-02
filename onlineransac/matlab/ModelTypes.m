function names = ModelTypes(withoutbern)
% Return a cell with as many model types as they exist, in the order they
% should be used for modelling.
%
% NAMES <- a cell with a text for each model

    names = { 'LL3', 'LN3', 'EXP2' };
    if ~withoutbern
        names{length(names)+1} = 'BERN';
    end

end