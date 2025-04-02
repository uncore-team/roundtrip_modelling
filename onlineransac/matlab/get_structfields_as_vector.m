function data = get_structfields_as_vector(c,fs)
% Given a cell C that contains a struct in each element, collect all values
% of the chain of fields FS (a cell; see access_deepest_field()) and form a
% vector with those values.

    n = length(c);
    data = nan(1,n);
    for f = 1:n
        data(f) = access_deepest_field(c{f},fs);
    end

end