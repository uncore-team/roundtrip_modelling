function defs = PerformanceDef(trace)
% Return the definition of all performance measures.
%
% TRACE -> 1 to show in console all performances.
%
% DEFS <- a cell with an entry per performance measure; that entry is a
%         struct:
%           .name <- name of the performance
%           .fields <- a cell with the sequence of field names to reach the
%                      performance into any performance object 
%                      (see PerformanceEmpty()).           

    emptyperf = PerformanceEmpty();
    fieldstree = all_deepest_fields(emptyperf);
    n = length(fieldstree);
    defs = cell(1,n);
    if trace
        fprintf('PERFORMANCE MEASURES (%d):\n',n);
    end
    for f = 1:n
        path = fieldstree{f};
        for g = 1:length(path)
            if g == 1
                na = path{g};
            else
                na = sprintf('%s:%s',na,path{g});
            end
        end
        if trace
            fprintf('\t%s\n',na);
        end
        defs{f} = struct('name',na,'fields',{path});
    end

end


