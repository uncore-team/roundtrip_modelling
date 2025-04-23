function r = ModelHasOffset(ty)
% Given the type of model, return 1 if that model has a parameter that is
% an offset.

    if strcmp(ty,'BERN')
        r = 0;
    else
        r = 1;
    end

end