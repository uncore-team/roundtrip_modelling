% test of some negative incorrect values in lognormal

clear;
close all;
clc;

params = [52412.4758766685  -7.36338669048975  0.988883157726085];

sample = [52412.4762071286 52412.4761548176 52412.4767982572 52412.4760871042 52412.4771032025 52412.4772688331 52412.4776455135 52412.4764980342 52412.4770460101 52412.4796478477 52412.4765448288 52412.4764084387 52412.4795355106 52412.4759297172 52412.4762803437 52412.4828894050 52412.4762744532 52412.4762409855 52412.4760652049 52412.4772655571];

[ok, params2(1), params2(2), params2(3)] = LognormalFit(sample);
if ~ok
    warning('ln3: ko');
else
    figure;
    [hfreqs,hxs] = hist(sample,50);
    bar(hxs,hfreqs);
    hold on;
    grid;
    xs = linspace(params(1),max(sample),100000);
    ys = LognormalPdf(xs,params(1),params(2),params(3));
    ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
    plot(xs,ys,'b-');
    xs = linspace(params2(1),max(sample),100000);
    yes = LognormalPdf(xs,params2(1),params2(2),params2(3));
    yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
    plot(xs,yes,'r--');
end
[reject, stat, thresh] = LognormalGof(sample, params2(1), params2(2), params2(3), 0)