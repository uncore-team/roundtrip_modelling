% Tabulate the threshold for the GoF experimentally for different sample
% sizes.

clear;
close all;
clc;

% model to tabulate and general experiment parms
mtrue = ModelCreate('LN3');
parmsunknown = 2; % 1 - tabulate for unknown parameters that are deduced from the sample; 2- same except the offset; 0 - tabulate for the true params that generate the sample
samplesizes = [20:10:500]; %20:10:500; % samples sizes to tabulate
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

    % deduce the threshold for getting an ALPHA level from the experimental
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
                results(f) = x0 + ((1 - alpha) - cdfsofar) / inca;
            else
                % Simpler estimate that agrees more with D'Agostino for the case of known parms: 
                results(f) = (hxs(g)+hxs(g-1))/2;
            end
            fprintf('\t\tThreshold: %f with cdf-1 = %f and cdf = %f (mean %f)\n',...
                    results(f),cdfsofar,newcdfsofar,(cdfsofar + newcdfsofar)/2);
            if traceinternal
                fprintf('\t');
                toc(t0)
            end            
            break;
        end
        cdfsofar = newcdfsofar;
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
plot(samplesizes,results,'.');
grid;
xlabel('sample size');
ylabel('threshold for gof');
subplot(1,2,2);
drawHisto(results,'histogram','threshold for gof');

%%  --- tries to fit the EXP2-parms-from-sample results

if ~strcmp(mtrue.type,'EXP2')
    error('Cannot run this section for non-EXP2 distrib');
end
% NOTES:
%   The curve tends to +Inf at 0 and has a more abrupt change in its
%   decreasing ratio than a/(X^b), ending almost horizontal at around y =
%   mean(results(29:end)) = 1.3966. For a size >= 300, this can be a
%   suitable threshold.
%   With the double exponential it fits ok in the Curve Fitter app:
%       Startpoint: [1.625;-0.0188;1.4653;-0.0002] (calculated
%       heuristically -"optimized"- by the app based on characteristics of
%       the data such as multiple decay speeds, initial values, preliminary 
%       one-exponential adjustment, derivatives, etc; actually Matlab does
%       not specify how)

% My code:

x0 = [1.625;-0.0188;1.4653;-0.0002];
% extended to 1000:
x0 = [1.11887690765281 -0.0100950663828462 1.37405116340726 -2.61317137849573e-05];
fcn1 = @(b,t) b(1).*exp(t.*b(2))+b(3).*exp(t.*b(4));
[parms, fval] = fminsearch (@ (b) norm(results - fcn1(b, samplesizes)), x0)
plot(samplesizes,results,'.');
grid;
hold on;
plot(samplesizes,fcn1(parms,samplesizes),'r.-');

% The Matlab automatic fitting code (produces the same result):

[xData, yData] = prepareCurveData( samplesizes, results );
% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% extended to 1000: opts.StartPoint = [1.11887690765281 -0.0100950663828462 1.37405116340726 -2.61317137849573e-05];
opts.StartPoint = [1.62051288611057 -0.0188232367797166 1.46532687005696 -0.000154533652141461];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'results vs. samplesizes', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'samplesizes', 'Interpreter', 'none' );
ylabel( 'results', 'Interpreter', 'none' );
grid on

%% --- tries to fit the LN3-parms-from-sample results

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
fv = horizfit - results(1); % final value (forced) for the 1st order response
indm = round(length(xs1)/2); % middle point taken to calc an initial guess
xm = xs1(indm);
ym = ys1(indm);
bstart = -log(-ym/fv+1)/xm % fit deterministically a 1st order with that middle point and final value
fcn1 = @(b,x) fv * b(1) * (1 - exp(-b(2) * x)); % forced to end at fv
[parms, fval] = fminsearch (@ (b) norm(ys1 - fcn1(b,xs1)), [1;bstart]) % fit a better 1st order, but it will not end at fv exactly (it cannot be forced because fv is at x = Inf, not at x = inds1(end))

plot(xs1 + samplesizes(1),ys1 + results(1),'o');
grid; 
hold on; 
plot(xs1(end) + samplesizes(1),fv + results(1),'m*');
plot(xs1 + samplesizes(1),fcn1(parms,xs1) + results(1),'r:.');

% first part alternative - 5th deg poly
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
