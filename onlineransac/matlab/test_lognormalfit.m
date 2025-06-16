% Testing LN3 fitting and its inner details

clear;
clc;
close all;

rng(54);

samplesizes = [20:10:200]; % [20:200, 210:10:2000, 2020:100:10000, 10000];
numtests = 100000;

nsamplesizes = length(samplesizes);
percdiffs = nan(1,nsamplesizes);
tvds = nan(nsamplesizes,numtests); % Total Variation Distance (TVD) between the true distribution and the fitted one (https://en.wikipedia.org/wiki/Total_variation_distance_of_probability_measures)
for ind = 1:nsamplesizes

    samplesize = samplesizes(ind);
    fprintf('\n\n=========== SAMPLESIZE %d\n',samplesize);

    casesigns = nan(1,numtests);
    stats = nan(1,numtests);
    numunfit = 0;
    oldperc = 0;
    t0 = tic;
    for t = 1:numtests
    
        perc = t / numtests * 100;
        if (round(perc/10) > round(oldperc/10))
            oldperc = perc;
            fprintf('\t%.2f%%, test %d out of %d. ',perc,t,numtests);
            toc(t0)
            drawnow;
        end
    
        finish = 0;
        while ~finish
            mo = ModelCreateRnd('LN3','typrnd'); % create the true params for the model (randomly)
            ds = ModelRnd(mo,1,samplesize); % draw a sample of the given size
            [ok, offs, mu, sigma, casesign] = LognormalFit(ds);
            if ok
                finish = 1;
            else
                numunfit = numunfit + 1;
            end
        end
        casesigns(t) = casesign;

        %[~,stat,~] = LognormalGof(ds,offs,mu,sigma,0);
        %stats(t) = stat;
        
        coeffstrue = ModelToCoeffs(mo);
        coeffstrue(2) = 0; % shift the distrib to 0 to avoid large numbers in the support for the trapz to be accurate
        mo0 = ModelFromCoeffs(coeffstrue);
        xs = linspace(0,ModelToExpectation(mo0,[]) + sqrt(ModelToVariance(mo0,[]))*2,10000);
        ysfit = LognormalPdf(xs,0,mu,sigma);
        ystrue = ModelPdf(mo0,[],xs);
        absysdiff = abs(ysfit-ystrue);
        tvds(ind,t) = trapz(xs,absysdiff)/2; % TVD = 0 → identical distributions; TVD = 1 → completely disjoint
        if (tvds(ind,t) < 0) || (tvds(ind,t) > 1)
            error('Invalid calculation of TVD');
        end
    end
    nsame = length(find(casesigns == 0));
    ndiff = length(find(casesigns == 1));
    percdiffs(ind) = ndiff/(nsame+ndiff);

    
    fprintf('Numtests: %d\n',numtests);
    fprintf('Numunfit:%d\n',numunfit);
    fprintf('Num same signs:%d (%.2f%%)\n',nsame,nsame/(nsame+ndiff));
    fprintf('Num diff signs:%d (%.2f%%)\n',ndiff,ndiff/(nsame+ndiff));
    
%    histogram(stats);
%    title('Histogram of GoF statistics');

end

figure;
plot(samplesizes,percdiffs,'.-');
title('Perc. of times there is a zero-cross in offset MLE');
xlabel('sample size');

figure;
boxplot(tvds.');
title('Abs. diff. between true and fitted distr.');
xlabel('sample size index');