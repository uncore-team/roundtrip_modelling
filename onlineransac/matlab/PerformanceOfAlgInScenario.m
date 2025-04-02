function res = PerformanceOfAlgInScenario(S,...
                                          algname,algfun,algparms,algstate0,...
                                          req,gtfun,perfparms,trace)
% Model the entire scenario S with the given algorithm and get performance 
% measures of that modelling.
%
% S -> vector of values (usually rtts) of the scenario.
% ALGNAME -> text indicating the algorithm to use.
% ALGFUN -> function of the algorithm (see AlgOnlineransac, AlgBernoulli, 
%           ...).
% ALGPARMS -> struct with the parameters of the algorithm.
% ALGSTATE0 -> struct with the initial state of the algorithm (see the
%              algorithm function).
% REQ -> requirement to satisfy by the system, or empty if none. See
%        ReqCreate.            
% GTFUN -> a function handler to get the ground-truth model for a given
%          index in the scenario. prototype: m = gtfun(index). Empty if no
%          gt model.
% PERFPARMS -> parameters for evaluating the performance of modelling (see
%              PerformanceCalc).
% TRACE -> 1 to trace in console.
%
% RES <- a cell with the results of modelling: one element per iteration
%        (empty if no iteration) where a model of the scenario was 
%        available. It is a struct:
%           .index -> index of the scenario where the model was available.
%           .perf -> performance struct returned by PerformanceCalc() for
%                    that index.

    n = length(S);
    res = {};
    if n <= 0
        return;
    end

    if trace
        fprintf('Modeling scenario of %d data with [%s]...\n',n,algname);
    end

    state = algstate0;
    oldperc = 0;
    t0 = tic;
    for f = 1:n
        if trace
            perc = f / n * 100;
            if (round(perc) > round(oldperc))
                oldperc = perc;
                fprintf('%.2f%%, f=%d, ',perc,f);
                toc(t0)
            end
        end
    
        t1 = tic;
        [~,state] = algfun(algparms,state,S(f));
        ct = toc(t1);
        if ~isempty(state.model)
            if ~isempty(gtfun)
                gtm = gtfun(f);
            else
                gtm = [];
            end
            res{length(res)+1} = struct('index',f,...
                                        'perf',PerformanceCalc(state.model,...
                                                               S,f,req,gtm,ct,perfparms));
        end
    end        
    
end