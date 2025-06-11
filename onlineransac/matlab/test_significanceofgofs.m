% test of the actual significance obtained by gof tests with and without
% fitting.

clc;
close all;
clear; 

% LL3
% mtrue = ModelCreate('LL3');
% mtrue = ModelChangeParms(mtrue,...
%                          [37086.059443685,3720.633937275,0.140459127]); % high prob of bad alpha
% bad alphas when too normal shape (c < 0.3)
%                         [0.1,5,0.05]); % high prob of bad alpha
%                         [0.1,5,0.25]); % very good estimate of alpha

% % LN3
% mtrue = ModelCreate('LN3');
% mtrue = ModelChangeParms(mtrue,...
%                          [40000,10,5]); % ignored if we create randomly the model

% EXP2
mtrue = ModelCreate('EXP2');
mtrue = ModelChangeParms(mtrue,...
                         [0,0]);


[numparms,namparms] = ModelParmsDef(mtrue.type);
numtests = 10000;
samplesize = 500; % max 500 for Gof thresholds have only be calculated up to that
fixedtrue = 0; % 0- use random models each time; 1- use always the true model
withfigs = 0;

suponiendoparms = 0; % # of rejects if we know true parms
nosuponiendoparms = 0; % # of rejects if we take parms from the sample
numunfit = 0;
historynosupparms = [];
stats = [];
if withfigs
    fi = figure;
end
rng(54);
for t = 1:numtests
       
    if ~fixedtrue
        mtrue = ModelCreateRnd(mtrue.type,'typrnd');
    end

    fprintf('%d... GT-model: ',t);
    ModelPrint(mtrue);

    ds = ModelRnd(mtrue,1,samplesize); % generate sample
    
    [reject1,stat1,~] = ModelGof(mtrue,ds,1); % gof with parms coming out of the sample
    fprintf('\trej:%d\n',reject1);
    suponiendoparms = suponiendoparms + reject1;

    mfit = ModelFit(ds,1,length(ds),mtrue.type);
    reject2 = 2;
    stat2 = NaN;
    if mfit.defined
        [reject2,stat2,thresh] = ModelGof(mfit,ds,0); % gof with parms coming from the sample
        fprintf('\tES-model:');
        ModelPrint(mfit);
        fprintf('\trej:%d\n',reject2);
        nosuponiendoparms = nosuponiendoparms + reject2;
        coeffstrue = ModelToCoeffs(mtrue);
        coeffsfit = ModelToCoeffs(mfit);
    else
        numunfit = numunfit + 1;
        fprintf('\tES-model: UNDEFINED\n');
    end
    historynosupparms = [historynosupparms ; ...
                         coeffstrue(2:2 + numparms - 1),coeffsfit(2:2 + numparms - 1),reject1,reject2];
    stats = [stats; ...
             stat1,stat2];
    
    if withfigs
        hold off;
        [hfreqs,hxs] = hist(ds,50);
        if reject1 == 1
            ec = 'r';
        else
            ec = 'c';
        end
        if reject2 == 1
            fc = 'y';
        else
            fc = 'g';
        end
        bar(hxs,hfreqs,'EdgeColor',ec,'FaceColor',fc);

        hold on;
        grid;
        xs = linspace(min(ds),max(ds),1000);
        ys = ModelPdf(mtrue,ds,xs);
        ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
        plot(xs,ys,'b-');

        if mfit.defined
            yes = ModelPdf(mfit,ds,xs);
            yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
            plot(xs,yes,'r--');    
        end

        title(sprintf('%d out of %d - rejects: %d  %d',t,numtests,reject1,reject2));
        drawnow;
    end
end

fprintf('\n\nALPHA (TYPE I ERROR) ESTIMATES:\n');

fprintf('Assuming parameters known:\n');
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/(numtests - numunfit));
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/(numtests - numunfit));
fprintf('\n');

fprintf('Assuming parameters unknown (undefined: %d; %.2f%%):\n',numunfit,numunfit/numtests*100);
fprintf('\tEst.alpha (Type I error): %f\n',nosuponiendoparms/(numtests - numunfit));
fprintf('\tCorrect detection: %f\n',1-nosuponiendoparms/(numtests - numunfit));


indrej1 = numparms * 2 + 1;
indrej2 = numparms * 2 + 2;
indscoeffstrue = 1:numparms;
indscoeffsfit = numparms + 1 : numparms * 2;

figure;
indsrej = find(historynosupparms(:,indrej2) == 1);
histrej = historynosupparms(indsrej,indscoeffsfit);
if numparms == 3
    plot3(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),histrej(:,indscoeffstrue(3)),'*r')
elseif numparms == 2
    plot(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),'*r')
else
    error('Invalid num of coeffs');
end
grid
hold on
indsnorej = find(historynosupparms(:,indrej2) == 0);
histnorej = historynosupparms(indsnorej,indscoeffsfit);
if numparms == 3
    plot3(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),histnorej(:,indscoeffstrue(3)),'.b')
else
    plot(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),'.b')
end
xlabel(namparms{1});
ylabel(namparms{2});
if numparms == 3
    zlabel(namparms{3});
end
title('Rejections with estimated fit');

if ~fixedtrue
    figure;
    grid
    hold on

    indsrej = find(historynosupparms(:,indrej1) == 1);
    histrej = historynosupparms(indsrej,indscoeffstrue);
    if numparms == 3
        plot3(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),histrej(:,indscoeffstrue(3)),'*r')
    else
        plot(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),'*r')
    end

    indsnorej = find(historynosupparms(:,indrej1) == 0);
    histnorej = historynosupparms(indsnorej,indscoeffstrue);
    if numparms == 3
        plot3(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),histnorej(:,indscoeffstrue(3)),'.b')
    else
        plot(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),'.b')
    end

    xlabel(namparms{1});
    ylabel(namparms{2});
    if numparms == 3
        zlabel(namparms{3});
    end
    title('Rejections with actual fit');
end

% Logistic Regression (Generalized Linear Model)
% Purpose: Estimate influence of each axis on the binary outcome.
% Interpretation:
%
%   Coefficients tell you the direction and magnitude of influence of each axis.
%   p-values assess the statistical significance of each predictor.
% Reference: Dobson & Barnett (2008), An Introduction to Generalized Linear Models.
fprintf('\n\nLOGISTIC REGRESSION (GENERALIZED LINEAR MODEL):\n');
varnames = namparms;
varnames{end + 1} = 'hyprej';
indsaux = find(historynosupparms(:,indrej2) ~= 2);
mdl = fitglm(historynosupparms(indsaux,indscoeffstrue),historynosupparms(indsaux,indrej2),...
             'Distribution', 'binomial',...
             'VarNames',varnames);
disp(mdl.Coefficients);


% Linear Discriminant Analysis (LDA)
% Purpose: Find a linear combination of features that separates the two classes.
% Interpretation:
%
%    The Linear field contains coefficients showing how each axis contributes to class separation.
% Reference: Fisher (1936), "The Use of Multiple Measurements in Taxonomic Problems", Annals of Eugenics.
fprintf('\n\nLINEAR DISCRIMINANT ANALYSIS:\n');
ldaModel = fitcdiscr(historynosupparms(indsaux,indscoeffstrue),historynosupparms(indsaux,indrej2));
disp(ldaModel.Coeffs(1,2).Linear)

% Mutual Information
% Purpose: Non-parametric measure of dependency between feature and class.
% MATLAB doesn't have mi built-in, but you can find implementations (e.g., FEX: mutualinfo by David Meyer).
% Reference: Cover & Thomas (2006), Elements of Information Theory.
%
% Requires Statistics and Machine Learning Toolbox
fprintf('\n\nMUTUAL INFORMATION:\n');
for f = 1:length(indscoeffstrue)
    mutualinfo = mi(historynosupparms(indsaux,indscoeffstrue(f)), historynosupparms(indsaux,indrej2));
    fprintf('Mutual info coeff #%d: %f\n',f,mutualinfo);
end


%%

% Deduce the threshold on the distribution of the statistics in order to
% get a 0.05 significance level of the test
indsvalid = find((~isinf(stats(:,2))) & (~isnan(stats(:,2))));
sts = stats(indsvalid,2);
[hfreqs,hxs] = hist(sts,100);
cdfsofar = 0.0;
a = trapz(hxs,hfreqs);
for f = 2:length(hxs)
    cdfsofar = cdfsofar + (hfreqs(f-1)+hfreqs(f))/a/2*(hxs(f)-hxs(f-1));
    if cdfsofar >= 1 - 0.05
        fprintf('Threshold: %f with cdf %f\n',(hxs(f)+hxs(f-1))/2,cdfsofar);
        break;
    end
end
figpdfs = figure;
hold on;
grid;
bar(hxs,hfreqs);
% non-parametrical spline fit (overfit!):
fitstatx = linspace(1e-9,max(sts),100000);
fitstaty = spline(hxs,hfreqs,fitstatx);
figure(figpdfs);
hspl = plot(fitstatx,fitstaty,'b.-');
% gamma fit:
gamfit = fitdist(sts,'gamma');
gampdf = pdf('gamma',fitstatx,gamfit.a,gamfit.b);
hgamma = plot(fitstatx,gampdf/trapz(fitstatx,gampdf)*trapz(hxs,hfreqs),'r.-');
% loglogistic fit:
mll3 = ModelFit(sts,1,length(sts),'LL3');
ModelPrint(mll3);
ll3fit = ModelPdf(mll3,sts,fitstatx);
figure(figpdfs);
hll3 = plot(fitstatx,ll3fit/trapz(fitstatx,ll3fit)*trapz(hxs,hfreqs),'m.-');
% lognormal fit:
mln3 = ModelFit(sts,1,length(sts),'LN3');
ModelPrint(mln3);
ln3fit = ModelPdf(mln3,sts,fitstatx);
figure(figpdfs);
hln3 = plot(fitstatx,ln3fit/trapz(fitstatx,ln3fit)*trapz(hxs,hfreqs),'c.-');
% exponential fit:
mex2 = ModelFit(sts,1,length(sts),'EXP2');
ModelPrint(mex2);
ex2fit = ModelPdf(mex2,sts,fitstatx);
figure(figpdfs);
hex2 = plot(fitstatx,ex2fit/trapz(fitstatx,ex2fit)*trapz(hxs,hfreqs),'g.-');
legend([hgamma,hspl,hll3,hln3,hex2],{'gamma/chi-sq','spline','LL3','LN3','EXP2'})
title('Distribution of the statistic (parms from the sample)');



% % fit a chi-square to the stats of the gof test with params coming from
% % data. A chi-square distribution with k degrees of freedom is equivalent to a 
% % gamma distribution with shape a=k/2 and scale b=2
% fprintf('Fitting of chi-square (gamma) to data:');
% pd = fitdist(stats(:,2), 'Gamma')
% thres = gaminv(1-0.05,pd.a,pd.b) % gives 43.653106943931576 for the LL3
% NOTE: EDF gof tests do not usually assume any standard distribution for
% the statistics; they tabulate as we do above


