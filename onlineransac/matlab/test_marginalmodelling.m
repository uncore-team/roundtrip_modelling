% Test of the scenarios characteristics when modelling with marginal
% distributions
% Each scenario is scanned with the three main distributions and the
% fitting (no GoF) results are gathered in order to have some idea of what
% parameters appear in our scenarios.

clear;
close all;
clc;

expcat = ExperimentCatalog(0);
numexps = length(expcat);

minsamplesize = 20;
maxsamplesize = 200;
stepsamplesize = 10;
discontinuoussegments = 1; % 0-> move segment window 1 datum each step; 1-> move segment window by its entire length, >1 -> move segment that length

numsamplesizes = maxsamplesize - minsamplesize + 1;
parmscollected = cell(3,... % EXP3, LN3, LL3
                      numsamplesizes);
                 % each element in the cell is a matrix of 3 x N or 2 x N 
                 % (depending on the number of parameters of the distrib), being
                 % N the number of segments scanned with that distribution
                 % in all scenarios
indexexp3 = 1;
indexln3 = 2;
indexll3 = 3;

totscans = 0;
measuredscans = [0;0;0]; % one per distribution
for indexp = 1:numexps

    expe = expcat{indexp};
    [~,~,data] = ExperimentGet(expcat,expe.class,expe.index,1,Inf,0,NaN,0);
    fprintf('\nExperiment %d out of %d: %s (length: %d)\n',indexp,numexps,expe.name,length(data));
    t0 = tic;
    for s = minsamplesize:stepsamplesize:maxsamplesize % different sample sizes

        fprintf('\ts = %d  ',s);
        toc
        if discontinuoussegments == 0
            st = 1;
        elseif discontinuoussegments == 1
            st = s;
        else
            if discontinuoussegments < 0 
                error('Invalid discontinuous segments')
            end
            st = discontinuoussegments;
        end
        for f = 1:st:(length(data)-s) % scanning window
            %fprintf('\t\twin[%d : %d]   ',f,f+s);
            %toc
            windata = data(f:f + s);

            [alpha,beta] = ExponentialFit(windata);
            parmscollected{indexexp3,s - minsamplesize + 1} = addparms(parmscollected{indexexp3,s - minsamplesize + 1}, [alpha;beta]);
            measuredscans(indexexp3) = measuredscans(indexexp3) + 1;

            [ok, offs, mu, sigma] = LognormalFit(windata);
            if ok
                parmscollected{indexln3,s - minsamplesize + 1} = addparms(parmscollected{indexln3,s - minsamplesize + 1}, [offs;mu;sigma]);
                measuredscans(indexln3) = measuredscans(indexln3) + 1;
            end

            [a, b, c, exitflag] = LoglogisticFit(windata);
            if exitflag >= 0
                parmscollected{indexll3,s - minsamplesize + 1} = addparms(parmscollected{indexll3,s - minsamplesize + 1}, [a;b;c]);
                measuredscans(indexll3) = measuredscans(indexll3) + 1;
            end

            totscans = totscans + 1;
        end

    end

end
save 'matlab_marginals.mat';



%% --- plotting results

close all;

showboxplot(indexexp3,1,'EXP2 - a','a',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);
showboxplot(indexexp3,2,'EXP2 - beta','beta',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);

showboxplot(indexln3,1,'LN3 - a','a',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);
showboxplot(indexln3,2,'LN3 - mu','mu',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);
showboxplot(indexln3,3,'LN3 - sigma','sigma',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);

showboxplot(indexll3,1,'LL3 - a','a',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);
showboxplot(indexll3,2,'LL3 - b','b',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);
showboxplot(indexll3,3,'LL3 - c','c',parmscollected,minsamplesize,stepsamplesize,maxsamplesize);


%% --- functions

function newparms = addparms(oldparms,parms)
% PARMS is column vector

    if isempty(oldparms)
        newparms = parms;
    else
        newparms = [ oldparms, parms ];
    end

end

function showboxplot(ind,row,tit,ylab,parmscollected,minsamplesize,stepsamplesize,maxsamplesize)

    boxplotdata = [];
    boxplotlabel = [];
    for s = minsamplesize:stepsamplesize:maxsamplesize % different sample sizes
        d = parmscollected{ind,s - minsamplesize + 1}(row,:).';
        boxplotdata = [boxplotdata; d];
        boxplotlabel = [boxplotlabel; repmat(s,length(d),1)];
    end
    figure;
    grid;
    hold on;
    boxplot(boxplotdata,boxplotlabel); %'PlotStyle','compact');
    title(sprintf('%s; min:%f,max:%f',tit,min(boxplotdata),max(boxplotdata)));
    xlabel('win.len.');
    ylabel(ylab);

end