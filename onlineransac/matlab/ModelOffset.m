function off = ModelOffset(mo)
% Given in MO a model, return its offset. If the model has no offset,
% return NaN.

    if ~mo.defined
        error('Cannot return offset of undefined model');
    end
    if ~ModelHasOffset(mo.type)
        off = NaN;
    else
        v = ModelToCoeffs(mo);
        off = v(2);
    end

end