% test of performance modelling in a scenario

clear;
close all;

% read the rtts from a file
data = load('rtts.txt');

% part of PerformanceOfAlgInScenario
algname = 'onlineransac'
%------

n = length(data);
if n <= 0
    return;
end

for modelpreserving = 0:2
    for samplesliding = 0:1
        for datapreserving = 0:1
            fprintf("model_preserving[%d] " + ...
                    "sample_sliding[%d] " + ...
                    "data_preserving[%d]", modelpreserving, samplesliding, datapreserving);
            onlineransacparms = struct('s',20,'N',NaN,...
                           'mtypes',{{'LL3'}},... % {{'LL3','LN3','EXP2'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
                           'modelpreserving',modelpreserving,...
                           'samplesliding',samplesliding,...
                           'datapreserving',datapreserving);

            state = struct('sample',[],'ind0',NaN,'ind1',NaN,'model',[],'consensus',[]);

            t0 = tic;
            for f = 1:n
                % t1 = tic;
                [~,state] = AlgOnlineransac(onlineransacparms,state,data(f));
                % ct = toc(t1);
            end
            fprintf("\nModel: %s\n" + ...
                    "\talpha: %f\n" + ...
                    "\tbeta: %f\n", state.model.type, state.model.coeffs.alpha, state.model.coeffs.beta);
        end
    end
end
