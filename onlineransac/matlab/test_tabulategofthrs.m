% Tabulate the threshold for the GoF experimentally for different sample
% sizes.

clear;
close all;
clc;

% model to tabulate and general experiment parms
mtrue = ModelCreate('LL3');
parmsunknown = 1; % 1 - tabulate for unknown parameters that are deduced from the sample; 0 - tabulate for the true params that generate the sample
samplesizes = 20:10:500; % samples sizes to tabulate
numtests = 10000; % monte carlo simulation on that number of samples
alpha = 0.05; % significance level to tabulate for

% internal parms that should not be changed
traceinternal = parmsunknown;
usesimplerinterp = 1; % 1 == adjust better to the thresholds in D'Agostino for known parms and produces more coherent results when tested with test_significanceofgof.m

measurements = cell(1,length(samplesizes));
results = nan(1,length(samplesizes));
t0ext = tic;
for f = 1:length(samplesizes)
    samplesize = samplesizes(f);
    fprintf('TABULATING %s for size %d... ',mtrue.type,samplesize);
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
            break;
        end
        cdfsofar = newcdfsofar;
    end
    
end

namefi = sprintf('matlab_tabulatethresgof_%s',mtrue.type);
if parmsunknown
    namefi = sprintf('%s_parmsfromsample',namefi);
else
    namefi = sprintf('%s_parmstrue',namefi);
end
save(sprintf('%s.mat',namefi));

figure;
subplot(1,2,1);
plot(samplesizes,results,'.');
grid;
xlabel('sample size');
ylabel('threshold for gof');
subplot(1,2,2);
drawHisto(results,'histogram','threshold for gof');