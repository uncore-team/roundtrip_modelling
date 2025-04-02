function n = get_uniques_in_cell(c)
% Given a cell C, return the another cell with no repetitions of elements.
%
% C -> a cell.
% 
% N <- cell with no repetitions.

    n = {};
    for f = 1:length(c)
        found = 0;
        for g = 1:length(n)
            if ischar(n{g})
                if strcmp(n{g},c{f}) 
                    found = 1;
                    break;
                end
            elseif (n{g} == c{f})
                found = 1;
                break;
            end
        end
        if ~found
            n{length(n) + 1} = c{f};
        end
    end

end