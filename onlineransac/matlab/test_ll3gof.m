% test of the gof of the LL3

clc;
close all;
clear; 

a = 0.1;
b = 5;
c = 0.25;
paramsyes = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.45]);
numtests = 100;
samplesize = 10000;

nosuponiendoparms = 0; % rejects if we take parms from the sample
suponiendoparms = 0; % rejects if we know true parms
falsosparms = 0;
fi = figure;
for t = 1:numtests
   
    fprintf('%d...\n',t);
    
    ds = LoglogisticRnd(paramsyes(1),paramsyes(2),paramsyes(3),1,samplesize);
    hold off;
    [hfreqs,hxs] = hist(ds,50);
    bar(hxs,hfreqs);
    hold on;
    grid;
    xs = linspace(0,max(ds)*1.1,1000);
    ys = LoglogisticPdf(xs,paramsyes(1),paramsyes(2),paramsyes(3));
    ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
    plot(xs,ys,'b-');

    [ae, be, ce, exitflag] = LoglogisticFit(ds);
    yes = LoglogisticPdf(xs,ae,be,ce);
    yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
    plot(xs,yes,'r--');
    [reject,stat,thresh] = LoglogisticGoF(ds,ae,be,ce,0);
    nosuponiendoparms = nosuponiendoparms + reject;
    
    [reject,stat,thresh] = LoglogisticGoF(ds,paramsyes(1),paramsyes(2),paramsyes(3),1);
    suponiendoparms = suponiendoparms + reject;

    paramsno = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.45]);
    [reject,stat,thresh] = LoglogisticGoF(ds,paramsno(1),paramsno(2),paramsno(3),1);
    falsosparms = falsosparms + reject;
    
    drawnow;
end

fprintf('Assuming parameters known:\n');
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/numtests);

fprintf('\n\n');

fprintf('Assuming parameters unknown:\n');
fprintf('\tEst.alpha (Type I error): %f\n',nosuponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-nosuponiendoparms/numtests);



%%

function num = rnd(xmin, xmax)
    num = xmin + (xmax - xmin) * rand();
end 

function params = generate_params(alim, blim, clim)

    if nargin == 0
        error('wrong params number!');
    end
    if nargin >= 1
        params(1) = rnd(alim(1), alim(2));
    end
    if nargin >= 2
        params(2) = rnd(blim(1), blim(2));
    end
    if nargin == 3
        params(3) = rnd(clim(1), clim(2));
    end
end 