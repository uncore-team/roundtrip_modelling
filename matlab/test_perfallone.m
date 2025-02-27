% Get one value of performance for every performance and every factor
% combination available

clear;
close all;

[fs,factornames,factorlevels,totlevs] = ExperimentFactorDefs(1);
numfactors = length(factornames);
perfdefs = PerformanceDef(1);
numperfs = length(perfdefs);
expcat = ExperimentCatalog(1);
numexps = length(expcat);
[algcat,nallcombs] = AlgCatalog(1);
numalgs = length(algcat);

% original data for Carmen: struct('prediction',struct('lookahead',20,'futuremodel','LL3'));
perfparms = struct('prediction',struct('lookahead',20,'futuremodel','LL3'));
indstartsc = 1;
indendsc = 200;

% scan all experiments with all algorithms and all factors
% collect a table with a row per index in each experiment that contains the
% levels of the factors used in the experiment plus the prediction expected
% error, using model for the future, at that index, i.e., the difference
% between the expectation of the current model and the expectation of the
% future model.

tableperfs = [];
fprintf('\n\nScanning all scenarios / all algorithms / all parameters...\n');
t0 = tic;
for indexp = 1:numexps
    expe = expcat{indexp};
    fprintf('Experiment %d out of %d: %s(%d,%d)\n',indexp,numexps,expe.name,indstartsc,indendsc);
    [~,~,data] = ExperimentGet(expe.class,expe.index,1,Inf,0,NaN,0);
    data = data(indstartsc:indendsc);
    meansc = mean(data);
    for indalg = 1:numalgs
        fprintf('    Algorithm %d out of %d: %s\n',indalg,numalgs,algcat{indalg}.name);
        for indpcs = 1:length(algcat{indalg}.parmscombs)
            fprintf('        Parms combination %d out of %d ', ...
                    indpcs,length(algcat{indalg}.parmscombs));
            toc(t0)

            perfs = PerformanceOfAlgInScenario(data,...
                                               algcat{indalg}.name,...
                                               algcat{indalg}.fun,...
                                               algcat{indalg}.parmscombs{indpcs},...
                                               algcat{indalg}.state0,...
                                               [],[],perfparms,0);

            if ~isempty(perfs)
                ind = length(perfs);
                while (ind > 1) && isnan(perfs{ind}.perf.prediction.experr.model)
                    ind = ind - 1;
                end
                if ind >= 1
                    [~,levsnum] = ExperimentFactors(expe.class,expe.index,...
                                                    algcat{indalg}.name,...
                                                    algcat{indalg}.parmscombs{indpcs});

                    % % just store the results at the last point of the
                    % % scenario where the performances were correctly
                    % % calculated, ignoring all previous points.
                    % tableperfs = [tableperfs ; ...
                    %               indexp,indalg,indpcs, ...                                  
                    %               levsnum, ...
                    %               perfs{ind}.perf.prediction.experr.model, ...
                    %               perfs{ind}.perf.prediction.experr.sample, ...
                    %               perfs{ind}.perf.prediction.experr.next, ...
                    %               perfs{ind}.perf.prediction.stderr.sample, ...
                    %               perfs{ind}.perf.prediction.stderr.model, ...
                    %               perfs{ind}.perf.prediction.modeerr.sample, ...
                    %               perfs{ind}.perf.prediction.modeerr.model, ...
                    %               perfs{ind}.perf.prediction.modeerr.next ];

                    % store the mean if all valid perf values for the entire scenario:
                    % divide errors by 'meansc' to normalize the error and
                    % make it relative: it makes sense to allow for more
                    % error if the roundtrip times are greater                    
                    meanperfs = zeros(1,8);
                    for in = 1:ind
                        meanperfs = meanperfs + [perfs{in}.perf.prediction.experr.model / meansc, ...
                                                  perfs{in}.perf.prediction.experr.sample / meansc, ...
                                                  perfs{in}.perf.prediction.experr.next / meansc, ...
                                                  perfs{in}.perf.prediction.stderr.sample, ...
                                                  perfs{in}.perf.prediction.stderr.model, ...
                                                  perfs{in}.perf.prediction.modeerr.sample / meansc, ...
                                                  perfs{in}.perf.prediction.modeerr.model / meansc, ...
                                                  perfs{in}.perf.prediction.modeerr.next / meansc];
                    end

                    meanperfs = meanperfs / ind;
                    tableperfs = [tableperfs ; ...
                                  indexp,indalg,indpcs, ...                                  
                                  levsnum, ...
                                  meanperfs ];

                    %fprintf('tableperfs: ind %d, %f\n',ind,perfs{ind}.perf.prediction.experr.model);
                end
            end
        end
    end

end