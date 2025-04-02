function perf = PerformanceCalc(m,S,index,req,gtm,ct,parms)
% Given a model M obtained at index INDEX of scenario S, calculate a number 
% of performance measures.
%
% M -> model (as in ModelCreate).
% S -> the complete scenario.
% INDEX -> the index in S where the model is currently assessed. It is
%          assumed that further indexes are in the future of the modelling
%          process.
% REQ -> time requeriment to satisfy or empty if none (see ReqCreate).
% GTM -> ground-truth model at current index, if available. Empty if none.
% CT -> computation time (secs) spent in modelling at this index.
% PARMS -> parameters for evaluation. A struct:
%           .prediction -> a struct with the prediction parameters:
%                             .lookahead -> number of indexes in the future 
%                                           to use for prediction 
%                                           performances.
%                             .futuremodel -> type of the model to fit in
%                                             the future sample (see
%                                             ModelFit).
%
% PERF <- a struct with the performance measures (see PerformanceEmpty).

    if isempty(m) || (~m.defined)
        error('Undefined model');
    end

    % calculate performances

    if ~isempty(req)
        probm = ReqCheckInModel(req,m,S,0);
    else
        probm = NaN;
    end
    if index < length(S)
        nextv = S(index + 1);
    else
        nextv = NaN;
    end
    expm = ModelToExpectation(m,S);
    stdm = sqrt(ModelToVariance(m,S));
    modem = ModelToMode(m,S);
    if index + parms.prediction.lookahead > length(S)
        futuresample = [];
        futuremodel = ModelCreate(parms.prediction.futuremodel);
    else
        ind0 = index + 1;
        ind1 = index + parms.prediction.lookahead;
        futuresample = S(ind0:ind1);
        futuremodel = ModelFit(futuresample,ind0,ind1,parms.prediction.futuremodel);
    end
    if ~isempty(futuresample)
        if ~isempty(req)
            probs = ReqCheckInSample(req,futuresample);
        else
            probs = NaN;
        end
        exps = mean(futuresample);
        stds = std(futuresample);
        modes = mode(futuresample);
    else
        probs = NaN;
        exps = NaN;
        stds = NaN;
        modes = NaN;
    end
    if futuremodel.defined
        if ~isempty(req)
            probfm = ReqCheckInModel(req,futuremodel,S,0);
        else
            probfm = NaN;
        end
        expfm = ModelToExpectation(futuremodel,S);
        stdfm = sqrt(ModelToVariance(futuremodel,S));
        modefm = ModelToMode(futuremodel,S);
    else
        probfm = NaN;
        expfm = NaN;
        stdfm = NaN;
        modefm = NaN;
    end
    if ~isempty(gtm)
        if ~isempty(req)
            probgt = ReqCheckInModel(req,gtm,S,0);
        else
            probgt = NaN;
        end
        expgt = ModelToExpectation(gtm,S);
        stdgt = sqrt(ModelToVariance(gtm,S));
        modegt = ModelToMode(gtm,S);
    end

    % fill performances

    perf = PerformanceEmpty();
    perf.computation.time = ct;    
    if ~isempty(gtm)
        perf.groundtruth.reqprob = probm - probgt;
        perf.groundtruth.experr = expm - expgt;
        perf.groundtruth.stderr = stdm - stdgt;
        perf.groundtruth.modeerr = modem - modegt;
    end
    if ~isempty(req)
        perf.prediction.reqerr.sample = probm - probs;
        perf.prediction.reqerr.model = probm - probfm;
    end
    perf.prediction.experr.sample = expm - exps;
    perf.prediction.experr.model = expm - expfm;
    %fprintf('%d %d %f %f %f\n',index,length(S),expm,expfm,perf.prediction.experr.model);
    perf.prediction.experr.next = expm - nextv;
    perf.prediction.stderr.sample = stdm - stds;
    perf.prediction.stderr.model = stdm - stdfm;
    perf.prediction.modeerr.sample = modem - modes;
    perf.prediction.modeerr.model = modem - modefm;
    perf.prediction.modeerr.next = modem - nextv;

end