% test of prediction behaviour
%
% -Utilidad de tener un modelo paramétrico en lugar del propio sample:
%   -Cuando se compara la performance de predicción, si el futuro se modela
%   en vez de usar la muestra de futuro concreta, el uso de un modelo del
%   pasado da muy buen resultado de predicción respecto a bernoulli. Cuando
%   el futuro no se modela para representar cualquier muestra que pudiera
%   haber venido, sino que se usa la muestra concreta, dan resultados
%   prácticamente iguales los modelos del pasado que el bernoulli del
%   pasado, ya que la muestra del futuro tiene el mismo comportamiento
%   estocástico que la del pasado, el cual se puede recoger bien con
%   Bernoulli; se pierde, de hecho la ventaja del siguiente punto, ya que
%   la muestra del futuro también tiene poca probabilidad de reflejar las
%   partes del soporte poco representadas en el pasado:
%   -Predicción de lo que pasa en los lugares del soporte donde no hay
%   sample pero sí modelo.
%   -Necesario el modelo para detectar cambios de régimen -> detecta
%   cuándo se producen los cambios de régimen.
%   -Posibilidad de reproducir el escenario aleatoriamente.

clear;
close all;

experimentbatch = 'realoct2023'; %'realoct2023'; %'sim'; % 'realoct2023' % 'realpapersensors'
experimentindex = 9;
futureismodel = 1; % 0-> future is raw sample; 1-> future is model from raw sample
experimentrange = [1,500];%length(historial)]; 
transformintorobotdistance = 0; % 1-> to transform rtts into the distance travelled by a robot that exhibits those rtts and moves with certain speed 
robotspeed = 0.25; % in m/s

[fdir,fexp,historial] = ExperimentGet(experimentbatch,experimentindex,...
                                      experimentrange(1),experimentrange(2),...
                                      transformintorobotdistance,robotspeed,...
                                      1);

if transformintorobotdistance
    units = 'meters';
else
    units = 'millisecs';
end
req = ReqCreate(0.9,...
                mean(historial) + std(historial),...
                mean(historial) + std(historial)*2);
thereisgt = 0;
if strcmp(experimentbatch,'sim') && (experimentindex == 1)
    % LoglogisticRnd(1000,2,0.3
    mgt = ModelCreate('LL3');
    mgt.coeffs.a = 1000;
    mgt.coeffs.b = 2;
    mgt.coeffs.c = 0.3;
    mgt.defined = 1;
    thereisgt = 1;
end

fprintf('Measuring prediction on %d round trip-times...\r\n',length(historial));
onlineransacparms = struct('s',20,'N',NaN,...
                           'mtypes',{{'LL3','LN3','EXP2'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
                           'modelpreserving',0,...
                           'samplesliding',1,...
                           'datapreserving',0);

n = length(historial);
lookahead = 20; % must be >= minimum needed to form models, if futureismodel
bernwin = 20;
ransacstate = struct('sample',[],'ind0',NaN,'ind1',NaN,'model',[],'consensus',[]);
models = {};
sizes = [];
result = []; % 1 -> index of check in the sample
             % 2 -> 1 if the model prob is >= the req probability             
             % 3 -> model prob
             % 4 -> bernoulli prob
             % 5 -> actual sample prob in the next LOOKAHEAD rtts
             % 6 -> underlying distrib. prob of satisfaction, if available
             % 7 -> model expectation
             % 8 -> bernoulli expectation
             % 9 -> actual sample expectation in the next LOOKAHEAD rtts
             % 10 -> model variance
             % 11 -> bernoulli variance
             % 12 -> actual sample variance in the next LOOKAHEAD rtts
             % 13 -> model mode
             % 14 -> bernoulli mode
             % 15 -> actual sample mode inthe next LOOKAHEAD rtts
oldperc = 0;
tic;
for f = 1:n-lookahead

    perc = f / n * 100;
    if (round(perc) > round(oldperc))
        oldperc = perc;
        fprintf('%.2f%%, f=%d, ',perc,f);
        toc
    end

    [exitbranch,ransacstate] = AlgOnlineransac(onlineransacparms,...
                                            ransacstate,...
                                            historial(f));
    if ~isempty(ransacstate.model)
        models{length(models) + 1} = struct('index',f,'model',ransacstate.model);
        sizes = [sizes length(ransacstate.sample)];

        pastsample = historial(max(1,f-bernwin-1):f);
        indfut0 = f + 1;
        indfut1 = f + lookahead;
        futuresample = historial(indfut0:indfut1);

        [okm,pm] = ReqCheckInModel(req,ransacstate.model,historial,0);        
        [~,psprev] = ReqCheckInSample(req,pastsample);
        if thereisgt
            [~,gtpr] = ReqCheckInModel(req,mgt,historial,0);
        else
            gtpr = NaN;
        end
        modelex = ModelToExpectation(ransacstate.model,historial);
        modelv = ModelToVariance(ransacstate.model,historial);
        modelmo = ModelToMode(ransacstate.model,historial);
        bernex = mean(pastsample);
        bernv = var(pastsample);
        bernmo = mode(pastsample);

        dowithmodel = futureismodel;
        if futureismodel            
            [m,~] = ModelAssess(futuresample,indfut0,indfut1,ModelTypes(1));
            if isempty(m) 
                dowithmodel = 0;
            else
                [~,ps] = ReqCheckInModel(req,m,historial,0);
                actex = ModelToExpectation(m,historial);
                actv = ModelToVariance(m,historial);
                actmo = ModelToMode(m,historial);
            end
        end
        if ~dowithmodel
            [~,ps] = ReqCheckInSample(req,futuresample);
            actex = mean(futuresample);
            actv = var(futuresample);
            actmo = mode(futuresample);
        end

        result = [result ; f,okm,...
                           pm,psprev,ps,gtpr,...
                           modelex,bernex,actex,...
                           modelv,bernv,actv,...
                           modelmo,bernmo,actmo];
    end
end

fprintf('Average time step: %f ms; std: %f ms\n',mean(historial),std(historial));

hfigesc = ScenarioShow(historial,fexp,[],models,units);
figure(hfigesc);
hold on;
ts = cumsum(historial);
plot([0 ts(end)],req.t0 * [1 1],'k--');
plot([0 ts(end)],req.t1 * [1 1],'k--');

figure;
plot(sizes,'.-');
grid;
title('Sample sizes of assessed samples');
ylabel('size');

if thereisgt
    showErrorComparison(historial,...
                   result(:,1),...
                   result(:,3)-result(:,6),...
                   result(:,4)-result(:,6),...
                   'Comparison to GT prob',...
                   'model-gt','bern-gt',...
                   units);
end
  
errsmodel = result(:,5)-result(:,3);
errsbern = result(:,5)-result(:,4);
showErrorComparison(historial,...
               result(:,1),...
               errsmodel,...
               errsbern,...
               'Comparison of distance to prob of req',...
               'model','bern',...
               units);
showErrorComparison(historial,...
               result(:,1),...
               result(:,7)-result(:,9),...
               result(:,8)-result(:,9),...
               'Comparison of expectations',...
               'E[model] - E[rtt]','E[bern] - E[rtt]',...
               units);
showErrorComparison(historial,...
               result(:,1),...
               result(:,10)-result(:,12),...
               result(:,11)-result(:,12),...
               'Comparison of variances',...
               'V[model] - V[rtt]','V[bern] - V[rtt]',...
               units);
showErrorComparison(historial,...
               result(:,1),...
               result(:,13)-result(:,15),...
               result(:,14)-result(:,15),...
               'Comparison of modes',...
               'Model[model] - Mode[rtt]','Mode[bern] - Mode[rtt]',...
               units);
showErrorComparison(historial,...
               result(:,1),...
               result(:,13)-historial(result(:,1)+1),...
               result(:,14)-historial(result(:,1)+1),...
               'Comparison mode-next-rtt',...
               'Model[model] - next-rtt','Mode[bern] - next-rtt',...
               units);


