% test of the gof of the LL3

clc;
close all;
clear; 

%paramsyes = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.45]);
paramsyes = [0.1,5,0.25]; % very good estimation of alpha
paramsyes = [0.1,5,0.05]; % high prob of bad estimation
paramsyes = [37086.059443685,3720.633937275,0.140459127]; % high prob of bad estimation
numtests = 1000;
samplesize = 10000;
fixedll3 = 0;
withfigs = 0;
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
       
    if ~fixedll3
        paramsyes = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.45]);
    end

    fprintf('%d...\tGT: (%.9f,%.9f,%.9f) ',t,paramsyes(1),paramsyes(2),paramsyes(3));

    ds = LoglogisticRnd(paramsyes(1),paramsyes(2),paramsyes(3),1,samplesize);
    if withfigs
        hold off;
        [hfreqs,hxs] = hist(ds,50);
        bar(hxs,hfreqs);
        hold on;
        grid;
    end
    
    [reject,stat,thresh] = LoglogisticGoF(ds,paramsyes(1),paramsyes(2),paramsyes(3),1);
    fprintf('  \trej:%d\n',reject);
    suponiendoparms = suponiendoparms + reject;
    if withfigs
        xs = linspace(0,max(ds)*1.1,1000);
        ys = LoglogisticPdf(xs,paramsyes(1),paramsyes(2),paramsyes(3));
        ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
        plot(xs,ys,'b-');
    end

    [ae, be, ce, exitflag] = LoglogisticFit(ds);
    if withfigs
        yes = LoglogisticPdf(xs,ae,be,ce);
        yes = yes / trapz(xs,yes) * trapz(hxs,hfreqs);
        plot(xs,yes,'r--');    
    end
    [reject,stat,thresh] = LoglogisticGoF(ds,ae,be,ce,0);
    fprintf('\tEst:(%.9f,%.9f,%.9f)  \trej:%d\n',ae,be,ce,reject);
    nosuponiendoparms = nosuponiendoparms + reject;

    historynosupparms = [historynosupparms ; paramsyes(1),paramsyes(2),paramsyes(3),ae,be,ce,reject];
    
    if withfigs
        drawnow;
    end
end

fprintf('Assuming parameters known:\n');
fprintf('\tEst.alpha (Type I error): %f\n',suponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-suponiendoparms/numtests);
fprintf('\n\n');

fprintf('Assuming parameters unknown:\n');
fprintf('\tEst.alpha (Type I error): %f\n',nosuponiendoparms/numtests);
fprintf('\tCorrect detection: %f\n',1-nosuponiendoparms/numtests);


figure;
indsrej = find(historynosupparms(:,7) == 1);
histrej = historynosupparms(indsrej,4:6);
plot3(histrej(:,1),histrej(:,2),histrej(:,3),'*r')
grid
hold on
indsnorej = find(historynosupparms(:,7) == 0);
histnorej = historynosupparms(indsnorej,4:6);
plot3(histnorej(:,1),histnorej(:,2),histnorej(:,3),'.b')
xlabel('a')
ylabel('b')
zlabel('c')
title('Rejections vs. estimated fit');

if ~fixedll3
    figure;
    indsrej = find(historynosupparms(:,7) == 1);
    histrej = historynosupparms(indsrej,1:3);
    plot3(histrej(:,1),histrej(:,2),histrej(:,3),'*r')
    grid
    hold on
    indsnorej = find(historynosupparms(:,7) == 0);
    histnorej = historynosupparms(indsnorej,1:3);
    plot3(histnorej(:,1),histnorej(:,2),histnorej(:,3),'.b')
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
mdl = fitglm(historynosupparms(:,1:3),historynosupparms(:,7),...
             'Distribution', 'binomial',...
             'VarNames',{'a','b','c','hyprej'});
disp(mdl.Coefficients);


% Linear Discriminant Analysis (LDA)
% Purpose: Find a linear combination of features that separates the two classes.
% Interpretation:
%
%    The Linear field contains coefficients showing how each axis contributes to class separation.
% Reference: Fisher (1936), "The Use of Multiple Measurements in Taxonomic Problems", Annals of Eugenics.
ldaModel = fitcdiscr(historynosupparms(:,1:3),historynosupparms(:,7));
disp(ldaModel.Coeffs(1,2).Linear)

% Mutual Information
% Purpose: Non-parametric measure of dependency between feature and class.
% MATLAB doesn't have mi built-in, but you can find implementations (e.g., FEX: mutualinfo by David Meyer).
% Reference: Cover & Thomas (2006), Elements of Information Theory.
%
% Requires Statistics and Machine Learning Toolbox
mutualinfo_A = mi(historynosupparms(:,1), historynosupparms(:,7)) 
mutualinfo_B = mi(historynosupparms(:,2), historynosupparms(:,7)) 
mutualinfo_C = mi(historynosupparms(:,3), historynosupparms(:,7)) 



%%

function num = rnd(xmin, xmax)
    num = xmin + (xmax - xmin) * rand();
end 

function params = generate_params(alim, blim, clim)

    if nargin == 0
        error('wrong params number!');
    end
    if nargin >= 1
        params(1) = rnd(alim(1), alim(2));
    end
    if nargin >= 2
        params(2) = rnd(blim(1), blim(2));
    end
    if nargin == 3
        params(3) = rnd(clim(1), clim(2));
    end
end 

