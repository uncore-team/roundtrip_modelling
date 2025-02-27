% get all scenarios data to save them

clear;
close all;

expcat = ExperimentCatalog(1);
numexps = length(expcat);

ds = cell(1,numexps);
for f = 1:numexps
    [~,~,data] = ExperimentGet(expcat{f}.class,expcat{f}.index,1,Inf,0,NaN,0);
    ds{f} = data;
end