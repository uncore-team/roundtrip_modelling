% test of performance modelling in a scenario

clear;
close all;

expclass = 'realoct2023';
expind = 12;
maxdatatouse = 800;%Inf;%800;
transformintorobotdistance = 0; % 1-> to transform rtts into the distance travelled by a robot that exhibits those rtts and moves with certain speed 
robotspeed = 0.25; % in m/s

% [fdir,fexp,historial] = ExperimentGet(expclass,expind,...
%                                       1,maxdatatouse,...
%                                       transformintorobotdistance,robotspeed,...
%                                       1);

% read the rtts from a file
load matlab_perfallone070724.mat

                           % 'mtypes',{{'LL3','LN3','EXP2'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
onlineransacparms = struct('s',20,'N',NaN,...
                           'mtypes',{{'LL3', 'EXP2', 'LN3'}},... % EXP2 tends to identify more, shorter, regimes, but worse ones; thus it should be used as secondary
                           'modelpreserving',0,...
                           'samplesliding',1,...
                           'datapreserving',0);
ransacstate = struct('sample',[],'ind0',NaN,'ind1',NaN,'model',[],'consensus',[]);
perfparms = struct('prediction',struct('lookahead',20,'futuremodel','LL3'));
% [levsstr,levsnum] = ExperimentFactors(expclass,expind,...
%                                       'onlineransac',onlineransacparms);

perfs = PerformanceOfAlgInScenario(data,...
                                   'onlineransac',@AlgOnlineransac,onlineransacparms,ransacstate,...
                                   [],...
                                   [],...
                                   perfparms,...
                                   1);

% fprintf('%d performance measurements gathered.\n',length(perfs));

% fprintf('Factor levels in this experiment:\n');
% levsstr
% levsnum

% vtcomp = PerformanceOne(perfs,{'perf','computation','time'});
% figure;
% drawHisto(vtcomp,'comptime','secs/iter');