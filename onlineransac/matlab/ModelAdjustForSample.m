function nm = ModelAdjustForSample(m,ds)
% Given a model M and a sample drawn from it in DS, return a new model of
% the same type that is adjusted for satisfying the requirements that may
% exist on the offset of the model (e.g., that no datum in the sample
% equals the offset).

global TOLROUNDTRIPS

    ConstantsInit();

    if ~m.defined
        error('Tried to adjust an undefined model');
    end

    nm = m;
    r = ModelHasOffset(m.type);
    if ~r
        return;
    end
    o = ModelOffset(m);

    mi = min(ds);
    if mi < o
        error('Invalid model for sample: min(sample) < offset');
    end
    switch r
        case 1 % model can produce samples with values == offset
            % nothing to adjust
        case 2 % model cannot produce such samples
            if mi == o 
                o = mi - TOLROUNDTRIPS/10; % separate the offset by one order of magnitude less than the tolerance
                nm = ModelChangeOffset(m,o);
            end
        otherwise
            error('Invalid kind of offset in model');
    end

end
