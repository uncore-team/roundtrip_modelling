function ind = ExperimentCatalogFind(ca,cl,in)
% Return the index in the catalog CA (produced by ExperimentCatalog) of the
% experiment of class CL and index IN within the class. Error if not found.

    for f = 1:length(ca)
        if strcmp(ca{f}.class,cl) && (ca{f}.index == in)
            ind = f;
            return;
        end
    end
    error('Experiment not found in catalog');

end