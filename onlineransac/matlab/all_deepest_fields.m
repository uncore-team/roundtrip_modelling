function result = all_deepest_fields(myStruct)
% Return a cell with one element per leaf of the struct structure; the
% element is a cell with one name per field in the path to that leaf.
    
    result = {};
    fields = fieldnames(myStruct);
    for f = 1:numel(fields)
        fieldname = fields{f};
        fieldvalue = myStruct.(fieldname);
        if isstruct(fieldvalue)
            deepres = all_deepest_fields(fieldvalue);
            for g = 1:length(deepres) % each one is a complete leaf until this field
                newlem = {};
                newlem{1} = fieldname;
                el = deepres{g};
                for h = 1:length(el)
                    newlem{length(newlem)+1} = el{h};
                end
                result{length(result)+1} = newlem;
            end
        else
            newlem = {};
            newlem{1} = fieldname;
            result{length(result)+1} = newlem; % add one element per leaf
        end
    end

end