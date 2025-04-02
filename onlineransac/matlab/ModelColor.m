function co = ModelColor(m)
% Return a different color for each model type.
%
% M -> model (see ModelCreate).
%
% CO <- model color (a letter)

    % red is reserved for erroneous models, for example
    
    if strcmp(m.type,'LL3')

        co = 'y';

    elseif strcmp(m.type,'EXP2')

        co = 'g';

    elseif strcmp(m.type,'LN3')

        co = 'c';

    elseif strcmp(m.type,'BERN')

        co = 'm';
        
    else
        error('Invalid model type');
    end    

end