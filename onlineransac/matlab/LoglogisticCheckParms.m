function LoglogisticCheckParms(a,b,c)
% Produce an error if the parms are invalid.

    if (a <= 0) || (c <= 0) || (c > 0.5) || (b <= 0)
        error('Invalid parameters for a loglogistic');
    end

end