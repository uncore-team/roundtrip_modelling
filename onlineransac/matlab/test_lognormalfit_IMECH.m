% Testing LN3 fitting and its inner details

clear;
clc;
close all;

rng(54);

% samplesizes = [20:10:2000, 2020:100:10000, 10000];
samplesizes = [20:200, 210:10:2000, 2020:100:10000, 10000];
numtests = 100000;

nsamplesizes = length(samplesizes);
percdiffs = nan(1,nsamplesizes);
for ind = 1:nsamplesizes

    samplesize = samplesizes(ind);
    fprintf('\n\n=========== SAMPLESIZE %d\n',samplesize);

    casesigns = nan(1,numtests);
    stats = nan(1,numtests);
    numunfit = 0;
    oldperc = 0;
    t0 = tic;
    parfor t = 1:numtests
    
        % perc = t / numtests * 100;
        % if (round(perc/10) > round(oldperc/10))
        %     oldperc = perc;
        %     fprintf('\t%.2f%%, test %d out of %d. ',perc,t,numtests);
        %     toc(t0)
        %     drawnow;
        % end
    
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

plot(samplesizes,percdiffs,'.-');

save test_lognormalfit.mat