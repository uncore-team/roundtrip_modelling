% Tabulate the threshold for the GoF experimentally for different sample
% sizes.

clear;
close all;
clc;

dists = {'LN3'}; %'EXP2', ,'LL3'};

% samplesizes = 1010:10:7000; %[20:10:500]; %510:10:1000; % samples sizes to tabulate
% samplesizes = [20:10:2000, 2020:100:10000, 10000]; % samples sizes to tabulate
% samplesizes = 20:200;
% samplesizes = 20:2000;
samplesizes = 680:1020;

numtests = 100000; % monte carlo simulation on that number of samples
alpha = 0.05; % significance level to tabulate for
traceinternal = 0;

ndists = length(dists);
nsizes = length(samplesizes);
blocksize = 100; % block size

for d = 1:ndists
dist = dists{d};

for parmsunknown = 0:1 % 1 - tabulate for unknown parameters that are deduced from the sample; 2- same except the offset; 0 - tabulate for the true params that generate the sample

% model to tabulate and general experiment parms
mtrue = ModelCreate(dist);

results = nan(1,nsizes);
measurements = nan(nsizes,numtests);

fileprefix = sprintf('tabulatethresgof_%s',dist);
if parmsunknown
    if parmsunknown == 2
        fileprefix = sprintf('%s_parmsfromsampleexceptoffset',fileprefix);
        nametrace = sprintf('%s (parms from sample except offset)',dist);
    else
        fileprefix = sprintf('%s_parmsfromsample',fileprefix);
        nametrace = sprintf('%s (parms from sample)',dist);
    end
else
    fileprefix = sprintf('%s_parmstrue',fileprefix);
    nametrace = sprintf('%s (parms true)',dist);
end

namefi = sprintf('%s_%d_to_%d.mat',fileprefix,samplesizes(1),samplesizes(end));

if exist(namefi,'file')
    continue;
end

x1 = 1;
x2 = blocksize;
state = 0;

for f = 1:nsizes
    samplesize = samplesizes(f);

    % check for saving data and being able of recovering it from them
     if (f >= x1) && (f <= x2) % iter 1: 1 to blocksize, ...
        if state == 0
            namefi_bak = sprintf('bak_%s_%d_to_%d_[%d-%d].mat',fileprefix,samplesizes(1),samplesizes(end), x1, x2);

            if exist(namefi_bak,'file')
                fprintf('\t\tRESTORING %s block [%d-%d] for size %d... ',nametrace,x1,x2,samplesize);

                data = load(namefi_bak,"stats_block","results_block");
                measurements(x1:x2,:) = data.stats_block;
                results(x1:x2) = data.results_block;

                fprintf('\trestored!\n');
                state = 1;
            else
                state = 2;
            end
        elseif state == 1 % skip this sample size
            continue;
        elseif state == 2
            % nothing to do here
        end

     else % backup data

        % was it previously restored?
        if state == 2
            stats_block = measurements(x1:x2,:);
            results_block = results(x1:x2);
            save(namefi_bak,"stats_block","results_block");
        end

        x1 = x1 + blocksize;
        x2 = x2 + blocksize;
        state = 0;
    end

    fprintf('TABULATING %s for size %d...\n',nametrace,samplesize);

    stats = nan(1,numtests);
    parfor t = 1:numtests

        finish = 0;
        while ~finish
            mo = ModelCreateRnd(dist,'typrnd'); % create the true params for the model (randomly)
            ds = ModelRnd(mo,1,samplesize); % draw a sample of the given size
            if parmsunknown
                mfit = ModelFit(ds,1,length(ds),dist); % fit a model to the sample
                if ~mfit.defined
                    continue;
                end
                if (parmsunknown == 2) && ModelHasOffset(dist)
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

    measurements(f,:) = stats;

%    results(f) = deducethresholdfromstatshist(stats,alpha,traceinternal);
    results(f) = deducethresholdfromstatsquantile(stats,alpha,traceinternal);

end

save(namefi, 'samplesizes', 'measurements', 'results');

end % parmsunknown
end % dists

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
