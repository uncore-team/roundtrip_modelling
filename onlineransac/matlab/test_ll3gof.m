% test of the gof of the LL3

clc;
close all;
clear; 

a = 0.1;
b = 5;
c = 0.25;
numtests = 100;
samplesize = 10000;

nosuponiendoparms = 0; % rejects if we take parms from the sample
suponiendoparms = 0; % rejects if we know true parms
fi = figure;
for t = 1:numtests
   
    fprintf('%d...\n',t);
    
    ds = LoglogisticRnd(a,b,c,1,samplesize);
    hold off;
    [hfreqs,hxs] = hist(ds,50);
    bar(hxs,hfreqs);
    hold on;
    grid;
    xs = linspace(0,max(ds)*1.1,1000);
    ys = LoglogisticPdf(xs,a,b,c);
    ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
    plot(xs,ys,'b-');

    [ae, be, ce, exitflag] = LoglogisticFit(ds);
    yes = LoglogisticPdf(xs,ae,be,ce);
    yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
    plot(xs,yes,'r--');
    [reject,stat,thresh] = LoglogisticGoF(ds,ae,be,ce,0);
    nosuponiendoparms = nosuponiendoparms + reject;
    
    [reject,stat,thresh] = LoglogisticGoF(ds,a,b,c,1);
    suponiendoparms = suponiendoparms + reject;
    
    drawnow;
end

fprintf('Assuming parameters known:\n');
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/numtests);

fprintf('\n\n');

fprintf('Assuming parameters unknown:\n');
fprintf('\tEst.alpha (Type I error): %f\n',nosuponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-nosuponiendoparms/numtests);