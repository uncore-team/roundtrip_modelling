function ind = find_string_in_cell(c,s)
% Find a string S into a cell of strings C. Return the index if found or
% NaN if not found. 
% If S is not a string, convert it to string to do the search.

    if ~ischar(s)
        s = num2str(s);
    end
    for f = 1:length(c)
        if strcmp(c{f},s)
            ind = f;
            return;
        end
    end
    ind = NaN;

end