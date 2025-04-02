function v = PerformanceOne(perfs,which)
% Given the result of PerformanceOfAlgInScenario() into PERFS and a cell
% with as many strings as needed to identify the deepest field of a given
% performance according to the perf definition in PerformanceDef(), return
% a row vector with as many numbers as gathered in PERFS for that
% performance measure.
% WHICH must begin with 'perf'.

    n = length(perfs);
    v = nan(1,n);
    for f = 1:n
        v(f) = access_deepest_field(perfs{f},which);
    end

end