% Show an scenario

clear;
close all;

transformintorobotdistance = 0; % 1-> to transform rtts into the distance travelled by a robot that exhibits those rtts and moves with certain speed 
robotspeed = 0.25; % in m/s

[fdir,fexp,historial] = ExperimentGet('realoct2023',12,...
                                      -Inf,Inf,...
                                      transformintorobotdistance,robotspeed, ...
                                      1);%'realoct2023',6,1);
fprintf('Data with %d roundtrip times.\n',length(historial));
if transformintorobotdistance
    units = 'meters';
else
    units = 'millisecs';
end

tit = sprintf('%s',fexp);
fprintf('Experiment figures...\r\n');
h = ScenarioShow(historial,tit,[],[],units);
figure(h);
hold on;
plot(cumsum(historial),movmean(historial,20),'m-','LineWidth',2);

fprintf('Done.\r\n');