function value = access_deepest_field(myStruct, fieldNames)

    value = myStruct;
    for i = 1:length(fieldNames)
        fieldName = fieldNames{i};
        value = value.(fieldName);
    end

end
