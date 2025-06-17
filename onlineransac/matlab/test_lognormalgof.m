% test of ln3 gof

clear;
clc;
close all;

mtrue = ModelCreate('LN3');
parmstrue = 1;
samplesizes = [20:100:10000];

rng(54);

figure;
for s = 1:length(samplesizes)

    fprintf('---> samplesize %d\n',samplesizes(s));

    numtests = 10000;
    stats = nan(1,numtests);
    numsans = zeros(1,numtests);
    oldperc = 0;
    t0 = tic;
    for f = 1:numtests
            
        perc = f / numtests * 100;
        if (round(perc/10) > round(oldperc/10))
            oldperc = perc;
            fprintf('\t%.2f%%, test %d out of %d. ',perc,f,numtests);
            toc(t0)
            drawnow;
        end

        mtrue = ModelCreateRnd(mtrue.type,'typrnd');
    
        gofalready = 0;
        while ~gofalready
            fitalready = 0;
            while ~fitalready
                ds = ModelRnd(mtrue,1,samplesizes(s));
                if ~parmstrue
                    mo = ModelFit(ds,1,length(ds),mtrue.type);
                    if mo.defined
                        fitalready = 1;
                    end
                else
                    mo = ModelAdjustForSample(mtrue,ds);
                    fitalready = 1;
                end
            end
            [reject,stat,thresh,numsanitized] = LognormalGof(ds,mo.coeffs.gamma,mo.coeffs.mu,mo.coeffs.sigma,parmstrue);
    %        [reject,stat,thres] = ModelGof(mo,ds,parmstrue); % parms true
            %if ~isnan(thresh) % otherwise, gof has not worked
                gofalready = 1;
            %end
        end
        stats(f) = stat;
        numsans(f) = numsanitized;
    
    end
    indsvalid = find(stats ~= inf);

    hold off;
    histogram(stats(indsvalid));
    hold on;
    plot(2.492 * [1 1],[0 1000],'r-');
    title(sprintf('SS: %d, parmstrue: %d;  quantile(0.95) = %f, invalid: %.2f%%, max.sanit.: %.2f%%, numtestssanit.: %.2f%%',...
                  samplesizes(s),...
                  parmstrue,...
                  quantile(stats(indsvalid),0.95),...
                  (numtests-length(indsvalid))/numtests*100,...
                  max(numsans)/samplesizes(s)*100,...
                  length(find(numsans > 0))/numtests * 100));
    drawnow;

end
