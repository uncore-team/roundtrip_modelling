% test of the actual significance obtained by gof tests with and without
% fitting.

clc;
close all;
clear; 

% LL3
mtrue = ModelCreate('LL3');
mtrue = ModelChangeParms(mtrue,...
                         [37086.059443685,3720.633937275,0.140459127]); % high prob of bad alpha
% bad alphas when too normal shape (c < 0.3)
%                         [0.1,5,0.05]); % high prob of bad alpha
%                         [0.1,5,0.25]); % very good estimate of alpha

% LN3
mtrue = ModelCreate('LN3');
mtrue = ModelChangeParms(mtrue,...
                         [40000,10,5]); % good alpha
% bad alphas when too normal shape
%                         [40000,-10,0.5]); % very bad alpha


% EXP2
mtrue = ModelCreate('EXP2');
mtrue = ModelChangeParms(mtrue,...
                         [0,0]);


[numparms,namparms] = ModelParmsDef(mtrue.type);
numtests = 1000;
samplesize = 10000;
fixedtrue = 0;
withfigs = 1;

suponiendoparms = 0; % rejects if we know true parms
nosuponiendoparms = 0; % rejects if we take parms from the sample
numunfit = 0;
historynosupparms = [];
if withfigs
    fi = figure;
end
for t = 1:numtests
       
    if ~fixedtrue
        mtrue = ModelCreateRnd(mtrue.type,'typrnd');
    end

    fprintf('%d... GT-model: ',t);
    ModelPrint(mtrue);

    ds = ModelRnd(mtrue,1,samplesize);
    
    [reject1,~,~] = ModelGof(mtrue,ds,1);
    fprintf('\trej:%d\n',reject1);
    suponiendoparms = suponiendoparms + reject1;

    mfit = ModelFit(ds,1,length(ds),mtrue.type);
    reject2 = 2;
    if mfit.defined
        [reject2,stat,thresh] = ModelGof(mfit,ds,0);
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
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/numtests);
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


