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

            onlineransacparms = struct('s',20,'N',NaN,...
                           'mtypes',{{'LN3'}},... % {{'LL3','LN3','EXP2'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
                           'modelpreserving',modelpreserving,...
                           'samplesliding',samplesliding,...
                           'datapreserving',datapreserving);

            state = struct('sample',[],'ind0',NaN,'ind1',NaN,'model',[],'consensus',[]);

            t0 = tic;
            for f = 1:n
                t1 = tic;
                [exitbranch,state] = AlgOnlineransac(onlineransacparms,state,data(f));
                ct = toc(t1);
                fprintf("#%d, exitbranch[%d], time[%f]\n", f, exitbranch, ct);
            end

            fprintf("model_preserving[%d] " + ...
                    "sample_sliding[%d] " + ...
                    "data_preserving[%d]", modelpreserving, samplesliding, datapreserving);
            if ~isempty(state.model) && state.model.defined
                switch state.model.type
                    case 'EXP2'
                        fprintf("\nModel: %s\n" + ...
                                "\talpha: %f\n" + ...
                                "\tbeta: %f\n", state.model.type, state.model.coeffs.alpha, state.model.coeffs.beta);
                    case 'LN3'
                        fprintf("\nModel: %s\n" + ...
                                "\tgamma: %f\n" + ...
                                "\tmu: %f\n" + ...
                                "\tsigma: %f\n", state.model.type, state.model.coeffs.gamma, state.model.coeffs.mu, state.model.coeffs.sigma);
                    otherwise
                        fprintf("\nUnknown error\n");
                end
            else
                fprintf("\nModel: NOT defined\n");
            end
        end
    end
end
