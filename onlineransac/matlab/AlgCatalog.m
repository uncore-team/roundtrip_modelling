function [algcat,nallcombs] = AlgCatalog(trace)
% Return the catalog of known algorithms to model scenarios.
%
% TRACE -> 1 to show in console the info about algorithms.
%
% ALGCAT <- a cell with as many entries as algorithms; each one a struct:
%           .name <- alg name.
%           .fun <- function handler for the algorithm.
%           .state0 <- initial state object for the function.
%           .parmscombs <- a cell with all the possibilities of algorithm
%                          parameters (those that do not vary are kept the
%                          same). Each element is a struct with the whole
%                          parameters object.
% NALLCOMBS <- total number of parms combinations for all algs.

    [fs,factornames,numfactorlevels,~] = ExperimentFactorDefs(0);
    nallcombs = 0;

    algcat = {};

    algcat{length(algcat)+1} = struct('name','onlineransac',...
                                      'fun',@AlgOnlineransac,...
                                      'state0',struct('sample',[],...
                                                      'ind0',NaN,'ind1',NaN,...
                                                      'model',[],...
                                                      'consensus',[]));
    parmsbase = struct('s',20,'N',Inf);
    algcat{end}.parmscombs = {};
    indmtypesasfactor = find_string_in_cell(factornames,'onlineransac-modelforms');
    indmodelpreservasfactor = find_string_in_cell(factornames,'onlineransac-modelpreserv');
    indsampleslidingasfactor = find_string_in_cell(factornames,'onlineransac-sampleslide');
    inddatapreservasfactor = find_string_in_cell(factornames,'onlineransac-datapreserv');
    for f1 = 1:numfactorlevels(indmtypesasfactor)
        mtypeslevel = fs{indmtypesasfactor}.levels{f1};
        if ~strcmp(mtypeslevel,'noransac')
            mtypesposs = strsplit(mtypeslevel,'-');
            parmsbase.mtypes = mtypesposs;
            for f2 = 1:numfactorlevels(indmodelpreservasfactor)
                if ~strcmp(fs{indmodelpreservasfactor}.levels{f2},'noransac')
                    parmsbase.modelpreserving = f2 - 2;
                    for f3 = 1:numfactorlevels(indsampleslidingasfactor)
                        if ~strcmp(fs{indsampleslidingasfactor}.levels{f3},'noransac')
                            parmsbase.samplesliding = f3 - 2;
                            for f4 = 1:numfactorlevels(inddatapreservasfactor)
                                if ~strcmp(fs{inddatapreservasfactor}.levels{f4},'noransac')
                                    parmsbase.datapreserving = f4 - 2;
                                    indalgcatcomb = length(algcat{end}.parmscombs)+1;
                                    algcat{end}.parmscombs{indalgcatcomb} = parmsbase;
                                end
                            end
                        end
                    end
                end
            end
        end
    end    
    nallcombs = nallcombs + length(algcat{end}.parmscombs);

    algcat{length(algcat)+1} = struct('name','bernoulli',...
                                      'fun',@AlgBernoulli,...
                                      'state0',struct('model',[],'ind0',NaN,'ind1',NaN));
    parmsbase = struct('s',20,'win',20);
    algcat{end}.parmscombs = {};
    algcat{end}.parmscombs{1} = parmsbase;
    nallcombs = nallcombs + length(algcat{end}.parmscombs);

    if trace
        fprintf('ALGORITHMS (%d):\n',length(algcat));
        for f = 1:length(algcat)
            fprintf('\t%s (%d parms combinations)\n',algcat{f}.name,length(algcat{f}.parmscombs));
        end
        fprintf('TOTAL PARMS COMBS: %d\n',nallcombs);
    end

end