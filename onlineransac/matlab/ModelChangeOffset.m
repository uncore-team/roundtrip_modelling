function nm = ModelChangeOffset(mo,newoff)
% If the model has offset, change it by NEWOFF without affecting the rest
% of its parameters.

    if ~mo.defined
        error('Cannot change offset in undefined model');
    end
    if ~ModelHasOffset(mo.type)
        error('Cannot change offset in no-offset model');
    end
    v = ModelToCoeffs(mo);
    v(2) = newoff;
    nm = ModelFromCoeffs(v);

end