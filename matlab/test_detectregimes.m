% test of detection of regimes

clear;
close all;

maxdatatouse = 800;
transformintorobotdistance = 0; % 1-> to transform rtts into the distance travelled by a robot that exhibits those rtts and moves with certain speed 
robotspeed = 0.25; % in m/s

[fdir,fexp,historial] = ExperimentGet('realoct2023',10,...
                                      1,maxdatatouse,...
                                      transformintorobotdistance,robotspeed,...
                                      1);%'realoct2023',6,1);
regs = [];
if length(historial) > maxdatatouse
    fprintf('Data truncated for regime detection.\n');
else
    fprintf('Data with %d roundtrip times.\n',length(historial));
end
if transformintorobotdistance
    units = 'meters';
else
    units = 'millisecs';
end

fprintf('Detecting regimes on %d round trip-times...\r\n',length(historial));
onlineransacparms = struct('s',20,'N',NaN,...
                           'mtypes',{{'LL3','LN3','EXP2'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
                           'modelpreserving',0,...
                           'samplesliding',1,...
                           'datapreserving',0);
regs = regimeDetection(historial,onlineransacparms,1);
fprintf('%d regimes detected.\n',size(regs,1));
sumtot = length(historial);
sumexpl = zeros(1,sumtot);
maxlenreg = 0;
for f = 1:size(regs,1)
    sumexpl(regs(f,1):regs(f,2)) = 1; % some regimes may overlap others (if sample sliding)
    l = regs(f,2)-regs(f,1)+1;
    if l > maxlenreg
        maxlenreg = l;
    end
end
sumexpl = sum(sumexpl);
fprintf('%d explained out of %d vs %d not explained (%.2f%% explained)\n',sumexpl,sumtot,sumtot-sumexpl,sumexpl/sumtot*100);
fprintf('Max reg len: %d\n',maxlenreg);
tit = sprintf('%s; modeltypes:%d; modelpreserv:%d; sliding:%d; datapreserv:%d',...
              fexp,...
              length(onlineransacparms.mtypes),...
              onlineransacparms.modelpreserving,...
              onlineransacparms.samplesliding,...
              onlineransacparms.datapreserving);

fprintf('Experiment figures...\r\n');
ScenarioShow(historial,tit,regs,[],units);

fprintf('Done.\r\n');