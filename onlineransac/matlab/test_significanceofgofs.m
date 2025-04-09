% test of the gof of the LL3

clc;
close all;
clear; 

%paramsyes = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.45]);
paramsyes = [0.1,5,0.25]; % very good estimation of alpha
paramsyes = [0.1,5,0.05]; % high prob of bad estimation
paramsyes = [37086.059443685,3720.633937275,0.140459127]; % high prob of bad estimation

mtrue = ModelFromCoeffs([0,...
                         0.1,5,0.25]);
numtests = 1000;
samplesize = 10000;
fixedtrue = 0;
withfigs = 1;
% Results with that fixed LL3 above:
%
% Assuming parameters known:
% 	Est.alpha (Type I error): 0.048400
% 	Correct detection: 0.951600
% 
% Assuming parameters unknown:
% 	Est.alpha (Type I error): 0.030300
% 	Correct detection: 0.969700
%
%
% Results if we use for each test a random LL3:
%
% Assuming parameters known:
% 	Est.alpha (Type I error): 0.046200
% 	Correct detection: 0.953800
%
% Assuming parameters unknown:
% 	Est.alpha (Type I error): 0.644500
% 	Correct detection: 0.355500
%
%
% Results in the point-per-point reject figure:
%   As c goes below 0.3, we get a lot of false rejections. This occurs in
%   very gaussian shapes.

suponiendoparms = 0; % rejects if we know true parms
nosuponiendoparms = 0; % rejects if we take parms from the sample
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
    if withfigs
        hold off;
        [hfreqs,hxs] = hist(ds,50);
        bar(hxs,hfreqs);
        hold on;
        grid;
    end
    
    [reject,~,~] = ModelGof(mtrue,ds,1);
    fprintf('\trej:%d\n',reject);
    suponiendoparms = suponiendoparms + reject;
    if withfigs
        xs = linspace(0,max(ds)*1.1,1000);
        ys = ModelPdf(mtrue,ds,xs);
        ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
        plot(xs,ys,'b-');
    end

    mfit = ModelFit(ds,1,length(ds),mtrue.type);
    if withfigs
        yes = ModelPdf(mfit,ds,xs);
        yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
        plot(xs,yes,'r--');    
    end
    [reject,stat,thresh] = ModelGof(mfit,ds,0);
    fprintf('\tES-model:');
    ModelPrint(mfit);
    fprintf('\trej:%d\n',reject);
    nosuponiendoparms = nosuponiendoparms + reject;
    coeffstrue = ModelToCoeffs(mtrue);
    coeffsfit = ModelToCoeffs(mfit);
    historynosupparms = [historynosupparms ; coeffstrue(2:end),coeffsfit(2:end),reject];
    
    if withfigs
        drawnow;
    end
end

fprintf('\n\nALPHA (TYPE I ERROR) ESTIMATES:\n');

fprintf('Assuming parameters known:\n');
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/numtests);
fprintf('\n');

fprintf('Assuming parameters unknown:\n');
fprintf('\tEst.alpha (Type I error): %f\n',nosuponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-nosuponiendoparms/numtests);


coeffstrue = ModelToCoeffs(mtrue);
numcoeffs = length(coeffstrue) - 1;
indrej = numcoeffs * 2 + 1;
indscoeffstrue = 1:numcoeffs;
indscoeffsfit = numcoeffs + 1 : numcoeffs * 2;

figure;
indsrej = find(historynosupparms(:,indrej) == 1);
histrej = historynosupparms(indsrej,indscoeffsfit);
if numcoeffs == 3
    plot3(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),histrej(:,indscoeffstrue(3)),'*r')
elseif numcoeffs == 2
    plot(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),'*r')
else
    error('Invalid num of coeffs');
end
grid
hold on
indsnorej = find(historynosupparms(:,indrej) == 0);
histnorej = historynosupparms(indsnorej,indscoeffsfit);
if numcoeffs == 3
    plot3(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),histnorej(:,indscoeffstrue(3)),'.b')
else
    plot(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),'.b')
end
xlabel('a')
ylabel('b')
zlabel('c')
title('Rejections vs. estimated fit');

if ~fixedtrue
    figure;
    indsrej = find(historynosupparms(:,indrej) == 1);
    histrej = historynosupparms(indsrej,indscoeffstrue);
    if numcoeffs == 3
        plot3(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),histrej(:,indscoeffstrue(3)),'*r')
    else
        plot(histrej(:,indscoeffstrue(1)),histrej(:,indscoeffstrue(2)),'*r')
    end
    grid
    hold on
    indsnorej = find(historynosupparms(:,indrej) == 0);
    histnorej = historynosupparms(indsnorej,indscoeffstrue);
    if numcoeffs == 3
        plot3(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),histnorej(:,indscoeffstrue(3)),'.b')
    else
        plot(histnorej(:,indscoeffstrue(1)),histnorej(:,indscoeffstrue(2)),'.b')
    end
    xlabel('a')
    ylabel('b')
    zlabel('c')
    title('Rejections vs. actual fit');
end

% Logistic Regression (Generalized Linear Model)
% Purpose: Estimate influence of each axis on the binary outcome.
% Interpretation:
%
%   Coefficients tell you the direction and magnitude of influence of each axis.
%   p-values assess the statistical significance of each predictor.
% Reference: Dobson & Barnett (2008), An Introduction to Generalized Linear Models.
fprintf('\n\nLOGISTIC REGRESSION (GENERALIZED LINEAR MODEL):\n');
mdl = fitglm(historynosupparms(:,indscoeffstrue),historynosupparms(:,indrej),...
             'Distribution', 'binomial',...
             'VarNames',{'a','b','c','hyprej'});
disp(mdl.Coefficients);


% Linear Discriminant Analysis (LDA)
% Purpose: Find a linear combination of features that separates the two classes.
% Interpretation:
%
%    The Linear field contains coefficients showing how each axis contributes to class separation.
% Reference: Fisher (1936), "The Use of Multiple Measurements in Taxonomic Problems", Annals of Eugenics.
fprintf('\n\nLINEAR DISCRIMINANT ANALYSIS:\n');
ldaModel = fitcdiscr(historynosupparms(:,indscoeffstrue),historynosupparms(:,indrej));
disp(ldaModel.Coeffs(1,2).Linear)

% Mutual Information
% Purpose: Non-parametric measure of dependency between feature and class.
% MATLAB doesn't have mi built-in, but you can find implementations (e.g., FEX: mutualinfo by David Meyer).
% Reference: Cover & Thomas (2006), Elements of Information Theory.
%
% Requires Statistics and Machine Learning Toolbox
fprintf('\n\nMUTUAL INFORMATION:\n');
for f = 1:length(indscoeffstrue)
    mutualinfo = mi(historynosupparms(:,indscoeffstrue(f)), historynosupparms(:,indrej));
    fprintf('Mutual info coeff #%d: %f\n',f,mutualinfo);
end


