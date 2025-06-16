% Tabulate the threshold for the GoF experimentally for different sample
% sizes.

% NOTE: To reduce the heat and consumption of the computer, if it has 24
% cores, being 16-23 of non-high performance, we can force in linux to use
% them only this way: taskset -c 16-23 /APLICACIONES/Matlab2023b/bin/matlab
% And them limit to 8 the number of cores to use in matlab with: maxNumCompThreads(8)

clear;
close all;
clc;

% model to tabulate and general experiment parms
mtrue = ModelCreate('EXP2');
parmsunknown = 1; % 1 - tabulate for unknown parameters that are deduced from the sample; 2- same except the offset; 0 - tabulate for the true params that generate the sample
samplesizes = [20:10:2000, 2020:100:10000, 10000]; %[20:10:500]; %510:10:1000; % samples sizes to tabulate
numtests = 100000; % monte carlo simulation on that number of samples
alpha = 0.05; % significance level to tabulate for

% internal parms that should not be changed
traceinternal = parmsunknown;
usesimplerinterp = 1; % 1 == adjust better to the thresholds in D'Agostino for known parms and produces more coherent results when tested with test_significanceofgof.m

measurements = cell(1,length(samplesizes));
results = nan(1,length(samplesizes));
if parmsunknown
    if parmsunknown == 2
        nametrace = sprintf('%s (parms from sample except offset)',mtrue.type);
    else
        nametrace = sprintf('%s (parms from sample)',mtrue.type);
    end
else
    nametrace = sprintf('%s (parms true)',mtrue.type);
end
t0ext = tic;
for f = 1:length(samplesizes)
    samplesize = samplesizes(f);
    fprintf('TABULATING %s for size %d... ',nametrace,samplesize);
    toc(t0ext)

    stats = nan(1,numtests);
    if traceinternal
        oldperc = 0;
        t0 = tic;
    end
    for t = 1:numtests
        if traceinternal
            perc = t / numtests * 100;
            if (round(perc/10) > round(oldperc/10))
                oldperc = perc;
                fprintf('\t%.2f%%, test %d out of %d. ',perc,t,numtests);
                toc(t0)
                drawnow;
            end
        end
        finish = 0;
        while ~finish
            mo = ModelCreateRnd(mtrue.type,'typrnd'); % create the true params for the model (randomly)
            ds = ModelRnd(mo,1,samplesize); % draw a sample of the given size
            if parmsunknown
                mfit = ModelFit(ds,1,length(ds),mtrue.type); % fit a model to the sample
                if ~mfit.defined
                    continue;
                end
                if (parmsunknown == 2) && ModelHasOffset(mtrue.type)
                    mfitcoeffs = ModelToCoeffs(mfit);
                    mtruecoeffs = ModelToCoeffs(mo);
                    mfitcoeffs(2) = mtruecoeffs(2);
                    mfit = ModelFromCoeffs(mfitcoeffs);
                end
                [~,stat,~] = ModelGof(mfit,ds,0); % get the statistic for that fitting
                if isinf(stat) || isnan(stat)
                    continue;
                end
            else
                [~,stat,~] = ModelGof(mo,ds,1); % get the statistic for the true parms
                if isinf(stat) || isnan(stat)
                    continue;
                end
            end
            stats(t) = stat;
            finish = 1;
        end
    end

    measurements{f} = stats;
    if traceinternal
        fprintf('\t\t');
    end
%    results(f) = deducethresholdfromstatshist(stats,alpha,traceinternal);
    results(f) = deducethresholdfromstatsquantile(stats,alpha,traceinternal);
    if traceinternal
        fprintf('\t');
        toc(t0)
    end            
    
end

namefi = sprintf('matlab_tabulatethresgof_%s',mtrue.type);
if parmsunknown
    if parmsunknown == 2
        namefi = sprintf('%s_parmsfromsampleexceptoffset',namefi);
    else
        namefi = sprintf('%s_parmsfromsample',namefi);
    end
else
    namefi = sprintf('%s_parmstrue',namefi);
end
save(sprintf('%s_%d_to_%d.mat',namefi,samplesizes(1),samplesizes(end)));

%%

figure;
subplot(1,2,1);
plot(samplesizes,results,'b.');
hold on;
ys = movmean(results,25);
plot(samplesizes,ys,'r.-');
grid;
xlabel('sample size');
ylabel('threshold for gof');
subplot(1,2,2);
drawHisto(results,'histogram','threshold for gof');

return;

%% --- modeling of the LN3 data (as in tabulatethresgof_LN3_parmsfromsample_20_to_10000_withzoomupto2000.mat)

% we pursue a good fitting with a small norm of residuals w.r.t. the ones
% of higher degrees and coefficients not too small in the polynomials
%
% Matlab fitting tool finds the minimization of squared errors, which in
% turn is the same as getting the expectation of the Y conditioned to the
% X, thus we are getting the expectation of the thresholds for each sample
% size, which is what we want.

close all;
clc;

part1end = find(samplesizes == 72,1);
coeffs1 = [-0.000737065644440   0.184544019197911   1.144649233557212]; % a quadratic fit done with the basic fitting of the figure; get a norm of residuals = 0.288
plotpolyfit(samplesizes,results,1,part1end,coeffs1,'part 1');

part2end = find(samplesizes == 120,1);
coeffs2 = [-0.000054441897049   0.018475635035727  -2.149691220487750  89.985450676534995]; % cubic fitting; normresid = 0.312693333622991
plotpolyfit(samplesizes,results,part1end,part2end,coeffs2,'part 2');

part3end = find(samplesizes == 680,1);
coeffs3 = [-0.000000000001180   0.000000002636083  -0.000002240327089   0.000918517509819  -0.174031641455638  14.974327904387987]; % 5th polynomial; normresid = 13.494
plotpolyfit(samplesizes,results,part2end,part3end,coeffs3,'part 3');

part4end = find(samplesizes == 1020,1); 
coeffs4 = [-0.0000007674684   0.0019230494045  -1.5058287142259   385.6552504121422]; % cubic fitting with normres = 97.6
plotpolyfit(samplesizes,results,part3end,part4end,coeffs4,'part 4');

part5end = length(samplesizes);
coeffs5 = [0.000000000044030  -0.000000977059096   0.051670442087110 -13.899828584965791]; % cubic fitting with normresid = 129.06
plotpolyfit(samplesizes,results,part4end,part5end,coeffs5,'part 5');

transweight = @(x,k,transitionpoint) 1 ./ (1 + exp(-k * (x - transitionpoint)));
ks = 0.05 * ones(1,numparts - 1); % transition weights from one part to the next
ks(1) = 1;
ks(2) = 0.5;

% --- blending of the polynomials

endsparts = [part1end,part2end,part3end,part4end,part5end];
coeffsparts = {coeffs1,coeffs2,coeffs3,coeffs4,coeffs5};
numparts = length(coeffsparts);

ninds = length(samplesizes);
x = samplesizes;
y = zeros(1,ninds);
w = zeros(1,ninds);
for f = 1:ninds

    findex = f;
    n = samplesizes(f); % this is the samplesize

    if findex <= part1end/2 % all are indexes on the vector samplesize
        parts = [1]; % polynomials to weld
    elseif findex <= (part1end + part2end)/2
        parts = [1,2];
        transitionpoint = samplesizes(part1end);
    elseif findex <= (part2end + part3end)/2
        parts = [2,3];
        transitionpoint = samplesizes(part2end);
    elseif findex <= (part3end + part4end)/2
        parts = [3,4];
        transitionpoint = samplesizes(part3end);
    elseif findex <= (part4end + part5end)/2
        parts = [4,5];
        transitionpoint = samplesizes(part4end);
    else
        parts = [5];
    end

    if length(parts) == 1
        y(f) = polyval(coeffsparts{parts(1)},x(f));
    else
        poly1 = polyval(coeffsparts{parts(1)},x(f));
        poly2 = polyval(coeffsparts{parts(2)},x(f));     
        w(f) = transweight(x(f),ks(parts(1)),transitionpoint);
        y(f) = poly1 * (1 - w(f)) + poly2 * w(f);
    end

end
figure;
grid;
hold on;
plot(samplesizes(1:endsparts(1)),results(1:endsparts(1)),'.');
for f = 2:numparts
    plot(samplesizes(endsparts(f-1):endsparts(f)),results(endsparts(f-1):endsparts(f)),'.');
    plot(ones(1,2) * samplesizes(endsparts(f-1)),[min(results),max(results)],'--');
end
plot(x,y,'o-');
%scatter(x,y,36,w,'filled');  % 36 is the marker size
%colormap(parula);
%colorbar;
xlabel('samplesizes');
ylabel('threshold');
title('blended fitting');



%%  --- OBSOLETE --- tries to model the EXP2-parms-from-sample in several ways

% NOTES:
%   The double exponential calculates the initial parms for the search
%       heuristically -"optimized"- by the matlab app based on characteristics of
%       the data such as multiple decay speeds, initial values, preliminary 
%       one-exponential adjustment, derivatives, etc; actually Matlab does
%       not specify how

close all;
modetofit = '2doubleexp'; % '1doubleexp', '2doubleexp'

if strcmp(modetofit,'1doubleexp')

    x0 = [1.625;-0.0188;1.4653;-0.0002]; % gathered from heuristics of the matlab fitting app
    fcn1 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
    [parms, fval] = fminsearch (@ (b) norm(results - fcn1(b, samplesizes)), x0)
    plot(samplesizes,results,'.');
    grid;
    hold on;
    ysfcn = fcn1(parms,samplesizes);
    plot(samplesizes,ysfcn,'r.-');
    
    figure;
    sss = (ysfcn - results).^2;
    histogram(sss,50);
    title(sprintf('residuals (squared errors) - mean %.15f',mean(sss)));

elseif strcmp(modetofit,'2doubleexp')
   
    funtomin = @(x) twodoubleexpsplittedformin(samplesizes,results,x);
    minx = samplesizes(floor(length(samplesizes)*0.01));
    maxx = samplesizes(ceil(length(samplesizes)*0.25));
    
    % find the optimum transition point between both exponentials by
    % minimizing the resulting mean square error between the exponentials
    % and the data

    % figure;
    % hold on;
    % grid;
    % for f = minx:maxx
    %     fval = funtomin(f);
    %     plot(f,fval,'r.');
    % end
    
    options = optimset('Display','iter');
    x = fminbnd(funtomin,...
                minx,...
                maxx,...
                options);
    transitionpoint = round(x);
    [sss,parms1,parms2,k] = twodoubleexpsplitted(samplesizes,results,transitionpoint);
    format long;
    fprintf('Transition point at samplesize = %d\n',transitionpoint);
    fprintf('Transition welding ct: %f\n',k);
    fprintf('parms1: \n');
    disp(parms1);
    fprintf('parms2: \n');
    disp(parms2);

    figure;
    grid;
    hold on;
    % first part
    inds1 = find(samplesizes <= transitionpoint);
    ss1 = samplesizes(inds1);
    rs1 = results(inds1);
    fcn1 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
    plot(ss1,rs1,'b.');
    xs1 = linspace(ss1(1),ss1(end),1000);
    ysfcn1 = fcn1(parms1,xs1);
    plot(xs1,ysfcn1,'r-','LineWidth',2);
    % second part
    inds2 = find(samplesizes > transitionpoint);
    ss2 = samplesizes(inds2);
    rs2 = results(inds2);
    fcn2 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
    plot(ss2,rs2,'c.');
    xs2 = linspace(ss2(1),ss2(end),1000);
    ysfcn2 = fcn2(parms2,xs2);
    plot(xs2,ysfcn2,'m-','LineWidth',2);
    % smoothed transition (continuous - C^inf) through a logistic function
    % of the weigths centered at the transition point
    transweight = @(x) 1 ./ (1 + exp(-k * (x - transitionpoint)));
    transfunc = @(b1,b2,x) fcn1(b1,x) .* (1 - transweight(x)) + fcn2(b2,x) .* transweight(x);
    xst = linspace(samplesizes(1),samplesizes(end),10000);
    yst = transfunc(parms1,parms2,xst);
    plot(xst,yst,'k--','LineWidth',3);

    figure;
    histogram(sss,50);
    title(sprintf('residuals (squared errors) - mean %.15f',mean(sss)));

else
    error('unimplemented exp threshold modelling for those samplesizes');
end



%% --- OBSOLETE --- tries to fit the LN3-parms-from-sample results

% NOTES:
%   No need for fit. Uncorrelated threshold approximately = 1.178 with
%   sigma = 0.186, i.e., 95% between 0.992 and 1.364

drawHisto(results,'LN3-parms-from-sample','MonteCarlo threshold');
mean(results)

%% --- tries to fit the LL3-parms-from-sample results

% NOTES:
%   Follow something similar to a log curve, but estabilizes for x > 100, 
%   with certain slight decay for x > 420.
%   It seems reasonable to divide the curve into these parts:
%       x < 120 -> log fit after some shift/scale; alternatively, fit a
%                  first order response; alternatively, a 5th deg poly. 
%                  Anyway, the final value should be around
%                  6.238 to link with the next part.
%       (a line interpolation between 110 and 120 to close the gap; visually,
%        a linear interpolation between 90 sand 120 is more acceptable)
%       x in [120,370] -> horizontal; mean value = 6.238, sigma = 0.120
%       x > 370 -> linear (decreasing) fit with fixed intercept to join the
%                  previous part at 6.238; i.e., y = ax+b where y = 6.238
%                  if x = 370 => b = 6.238 - a * 370
%                  This is done this way:
%                       - substract 6.238 from the data to have an intercept
%                       of 0 in the new data y'
%                       - least-squares fit a curve y = mx to the new data by
%                       this: m = sum(x .* y') / sum(x.^2);
%                       - the resulting fit is y = mx + 6.238

close all;
clc;
figure;
grid;
hold on;

% second part

inds2 = find((samplesizes >= 120) & (samplesizes <= 370));
xs2 = samplesizes(inds2);
ys2 = results(inds2); 
horizfit = mean(ys2) % horizontal fit
ys2fit = horizfit * ones(1,length(xs2));
plot(xs2,ys2,'mo');
plot(xs2,ys2fit,'r.-');

% first part 

inds1 = find(samplesizes < 120);
xs1 = samplesizes(inds1) - samplesizes(1); % shift data to make x_0 = 0
ys1 = results(inds1) - results(1); % shift data to make y(0) = 0
plot(xs1 + samplesizes(1),ys1 + results(1),'o');
grid; 
hold on; 
% % fitting of a 1st order
% fv = horizfit - results(1); % final value (forced) for the 1st order response
% indm = round(length(xs1)/2); % middle point taken to calc an initial guess
% xm = xs1(indm);
% ym = ys1(indm);
% bstart = -log(-ym/fv+1)/xm % fit deterministically a 1st order with that middle point and final value
% fcn1 = @(b,x) fv * b(1) * (1 - exp(-b(2) * x)); % forced to end at fv
% [parms, fval] = fminsearch (@ (b) norm(ys1 - fcn1(b,xs1)), [1;bstart]) % fit a better 1st order, but it will not end at fv exactly (it cannot be forced because fv is at x = Inf, not at x = inds1(end))
% plot(xs1(end) + samplesizes(1),fv + results(1),'m*');
% plot(xs1 + samplesizes(1),fcn1(parms,xs1) + results(1),'r:.');

% fitting a 5th deg poly
%
% f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
% 
% Coefficients and 95% Confidence Bounds
%        Value      Lower     Upper  
% p1    -0.0000    -0.0000    0.0000
% p2    0.0000     -0.0000    0.0000
% p3    -0.0000    -0.0002    0.0002
% p4    -0.0015    -0.0139    0.0109
% p5    0.1927     -0.1347    0.5201
% p6    -0.7079    -3.8399    2.4242
% 
% Goodness of Fit
%             Value  
% SSE         0.0117
% R-square    0.9992
% DFE         4.0000
% Adj R-sq    0.9982
% RMSE        0.0541

% p1 is x^5. It looks like it is rather a third or second order poly, but
% actually quitting those coeffs it does not work; fitting lesser polys
% works worst too
parms = [-0.000000000421424;0.000000118010409;-0.000004990075374;-0.001491034210112;0.192700180267995;-0.707864114472733];
xs1 = samplesizes(inds1);
ys1fit5 = parms(1) * xs1.^5 + parms(2) * xs1.^4 + parms(3) * xs1.^3 + parms(4) * xs1.^2 + parms(5) * xs1 + parms(6);
plot(xs1,ys1fit5,'r.-');

% automatically generated code from Curve Fitting app:
%
% [xData, yData] = prepareCurveData( xs1, ys1 ); % xs1,ys1 without shifting
% or scaling!
% 
% % Set up fittype and options.
% ft = fittype( 'poly5' );
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'ys1 vs. xs1', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'xs1', 'Interpreter', 'none' );
% ylabel( 'ys1', 'Interpreter', 'none' );
% grid on


% third part

inds3 = find(samplesizes > 370);
xs3 = samplesizes(inds3) - samplesizes(inds3(1));
ys3 = results(inds3) - horizfit;
m = sum(xs3 .* ys3) / sum(xs3.^2)
plot(xs3 + samplesizes(inds3(1)),ys3 + horizfit,'co');
plot(xs3 + samplesizes(inds3(1)),m * xs3 + horizfit,'r.-');

% ALL PARTS TOGETHER - code to go to LoglogisticGoF.m

xs = samplesizes(1):samplesizes(end);
ys = [];
for n = xs

    if n < 90 % firs part: 5th poly

        parms = [-0.000000000421424;0.000000118010409;-0.000004990075374;-0.001491034210112;0.192700180267995;-0.707864114472733];
        thresh = parms(1) * n^5 + parms(2) * n^4 + parms(3) * n^3 + parms(4) * n^2 + parms(5) * n + parms(6);

    elseif n < 120 % link between first and second part: linear interpolation

        x1 = 90;
        y1 = 6.174206416983618; % value of the first part at 90
        x2 = 120;
        y2 = 6.237947960378145; % value of the second part at 120
        thresh = (y2 - y1)/(x2 - x1) * (n - x1) + y1;

    elseif n < 380 % second part: horizontal

        thresh = 6.237947960378145;

    else % third part: decaying linear

        m = -0.002007806860578; % slope
        thresh = m * (n - 380) + 6.237947960378145; % line that starts at n = 380 (second part)

    end

    ys = [ys,thresh];
end
plot(xs,ys,'g-','LineWidth',2);


%% AUXILIARY FUNCTIONS

function threshold = deducethresholdfromstatsquantile(stats,alpha,traceinternal)
% Deduce the threshold counting the stats from smallest to largest (more
% robust than from largest to smallest), until reaching 1-alpha area.
% This method is called "Empirical Quantile estimation", being the quantile
% 0.95 in the case that alpha == 0.05.

    threshold = quantile(stats,1-alpha);

    % % theoretically equivalent to the following, but matlab also
    % % interpolates to reach the eact point:
    % 
    % ntotstats = length(stats);
    % if ntotstats <= 1
    %     error('Cannot find threshold with <= 1 stats');
    % end
    % 
    % % number of stats needed to represent 1-alpha proportion:
    % countforminusalpha = round((1 - alpha) * ntotstats);
    % if (countforminusalpha <= 0) || (countforminusalpha >= ntotstats)
    %     error('Cannot find threshold with invalid countforminusalpha');
    % end
    % 
    % sortstats = sort(stats);
    % threshold = sortstats(countforminusalpha + 1); % first stat that goes to alpha proportion
    % if traceinternal
    %     fprintf('Threshold: %f; countminusalpha: %d; total: %d; indexfirstalpha: %d\n',...
    %             threshold,countforminusalpha,ntotstats,countforminusalpha + 1);
    % end
    
    if traceinternal
        fprintf('Threshold: %f\n',threshold);
    end
    
end

function threshold = deducethresholdfromstatshist(stats,alpha,traceinternal)
% Deduce the threshold for getting an ALPHA level from the experimental
% distribution of the statistic. Use the left side of the distribution
% since it has more data and therefore it is more robust than the right
% tail.

    [hfreqs,hxs] = hist(stats,100); % HXS are the equidistant centers of data
    cdfsofar = 0.0;
    a = trapz(hxs,hfreqs);
    for g = 2:length(hxs)
        newcdfsofar = cdfsofar + (hfreqs(g-1)+hfreqs(g))/a/2*(hxs(g)-hxs(g-1));
        if newcdfsofar >= 1 - alpha
            if ~usesimplerinterp
                % INCX is (HXS(g) - HXS(g-1)), the inter-center distance
                incx = hxs(g) - hxs(g-1);
                % INCA is NEWCDFSOFAR - CDFSOFAR            
                inca = newcdfsofar - cdfsofar;
                % At x0 = HXS(g - 1) + INCX/2 there is CDFSOFAR (< 1 - alpha)
                x0 = hxs(g - 1) + incx / 2;
                % At x1 = HXS(g) + INCX/2 there is NEWCDFSOFAR (>= 1 - alpha)
                % Between x0 and x1 we assume a uniform distribution of the
                %  data, thus the area between them is a rectangle: INCX * INCA
                % Therefore, 1 - alpha is reached when a subarea within x0 & x1
                %  gets exactly (1 - alpha) - CDFSOFAR:
                %   d * INCA = (1 - alpha) - CDFSOFAR, with d in [0,INCX]
                %   d = ((1 - alpha) - CDFSOFAR) / INCA
                %   The threshold should be at x0 + d.
                threshold = x0 + ((1 - alpha) - cdfsofar) / inca;
            else
                % Simpler estimate that agrees more with D'Agostino for the case of known parms: 
                threshold = (hxs(g)+hxs(g-1))/2;
            end
            break;
        end
        cdfsofar = newcdfsofar;
    end
    if traceinternal
        fprintf('Threshold: %f with cdf-1 = %f and cdf = %f (mean %f)\n',...
                threshold,cdfsofar,newcdfsofar,(cdfsofar + newcdfsofar)/2);
    end

end

function [sss,parms1,parms2,k] = twodoubleexpsplitted(samplesizes,results,tpoint)

    % here, samplesizes and results are the concatenation of two
    % experiments, in the original tests, one from 20 to 500 and one from 510 to 1000. I have
    % done the concatenation manually in the console using the .mat files of
    % those experiments.

    % two different sections are visually distinguished; both can be fitted
    % with double exponentials. The division in the original tests was at samplesize == 200
    % if only one section covering everything is desired, the double
    % exponential can be fit by starting at x0 = [1.11887690765281 -0.0100950663828462 1.37405116340726 -2.61317137849573e-05];

    transitionpoint = round(tpoint);
    inds1 = find(samplesizes <= transitionpoint);
    inds2 = find(samplesizes > transitionpoint);
    ss1 = samplesizes(inds1);
    rs1 = results(inds1);
    ss2 = samplesizes(inds2);
    rs2 = results(inds2);
    % fitting first part
    x01 = [2.78176575411404;-0.0418139301125338;1.68968288428376;-0.00074419627955785];
    fcn1 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
    [parms1, ~] = fminsearch (@ (b) norm(rs1 - fcn1(b, ss1)), x01);
    % fitting second part
    x02 = [0.916029848965478;-0.0122508889315835;1.41972515320634;-5.26936347849917e-05];
    fcn2 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
    [parms2, ~] = fminsearch (@ (b) norm(rs2 - fcn2(b, ss2)), x02);
    % smoothed transition (continuous - C^inf) through a logistic function
    % of the weigths centered at the transition point
    k = 1e-1; % a parameter for the logistic to bend both curves only locally
    transweight = @(x) 1 ./ (1 + exp(-k * (x - transitionpoint)));
    transfunc = @(b1,b2,x) fcn1(b1,x) .* (1 - transweight(x)) + fcn2(b2,x) .* transweight(x);
    ysfcn = transfunc(parms1,parms2,samplesizes);
    sss = (ysfcn - results).^2;

end

function meansss = twodoubleexpsplittedformin(samplesizes,results,transitionpoint)
    
    [sss,~,~,~] = twodoubleexpsplitted(samplesizes,results,transitionpoint);
    meansss = mean(sss);

end

function plotpolyfit(samplesizes,results,index1,index2,polycoeffs,tit)

    figure;
    grid;
    hold on;
    title(sprintf('Fitting of thresholds for %s',tit));
    plot(samplesizes(index1:index2),results(index1:index2),'b.');
    ys = polyval(polycoeffs,samplesizes(index1:index2));
    plot(samplesizes(index1:index2),ys,'r-');
    xlabel('Sample sizes');
    ylabel('thresholds');
    legend('data','polyfit');

end
