function data = get_structfields_as_cell(c,fs,tostr)
% Given a cell C that contains a struct in each element, collect all values
% of the chain of fields FS (a cell; see access_deepest_field()) and form a
% cell with those values.
%
% TOSTR -> 1 to convert the values to string before adding them to the
%          resulting cell in case they are not strings; 0 to do nothing.

    n = length(c);
    data = cell(1,n);
    for f = 1:n
        v = access_deepest_field(c{f},fs);
        if tostr && (~ischar(v))
            data{f} = num2str(v);
        else
            data{f} = v;
        end
    end

end