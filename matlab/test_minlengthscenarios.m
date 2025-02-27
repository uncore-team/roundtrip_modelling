% get the minimum length for scenarios

clear;
close all;

ca = ExperimentCatalog(0);
n = length(ca);
mile = Inf;
for f = 1:n
    fprintf('Loading experiment [%s]...\n',ca{f}.name);
    [fdir,fexp,data] = ExperimentGet(ca{f}.class,ca{f}.index,...
                                     -Inf,Inf,...
                                     0,NaN,...
                                     0);
    l = length(data);
    if l < mile
        mile = l;
    end
end
fprintf('Minimum length scenario: %d\n',mile);