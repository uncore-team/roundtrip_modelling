% test for deducing which factors have stronger influence on the effect

clear;
close all;
clc;

load('matlab_perfallone270124.mat');

addpath('GPfun_matlab_port');

numfactorstoconsider = 10; % 10 for all factors

X = tableperfs(:,4:4 + numfactorstoconsider - 1);
y = tableperfs(:,14);
ns = find(isnan(y));
nns = setdiff([1:length(y)],ns);
X = X(nns,:);
y = y(nns,:);
n = length(y);
d = size(X,2);
[fs,factornames,factorlevels,totlevs] = ExperimentFactorDefs(0);


% Perform regression 
var_names = factornames(1:numfactorstoconsider);
tbl = array2table([X, y], 'VariableNames', [var_names, 'Response']);
mdl = fitlm(tbl); % Adjust the formula as per your reduced dimensions

% Display the regression model summary
disp(mdl);

% Order most relevant factors
[sortvals, idx] = sort(abs(mdl.Coefficients.Estimate(2:end)), 'descend');
fprintf('Most relevant dimensions:\n');
for f = 1:d
    ind = idx(f);
    fprintf(['\t#%d: %s\t, factor-weight (abs): %f,\tp-value: %f\n'],...
            f,factornames{ind},sortvals(f),mdl.Coefficients.pValue(ind + 1));
end

% ---- GP (mirar el completo que tengo en GPfun_matlab_port hecho con
% Picasso)

fprintf('WITH GP...\n');

data = tableperfs;
indsnan = find(isnan(data(:,14)));
indsnnan = setdiff([1:size(data,1)],indsnan);
data = data(indsnnan,:);

ndata = size(data,1);
nfactorstoconsider = 5;

X_train = data(:,4:4+nfactorstoconsider-1);
[C,IA,IC] = unique(X_train,'rows');
Y_train = zeros(size(C,1),1);
for f = 1:size(C,1)
    indsindata = find(IC == f);
    correspondingRows = data(indsindata,14);
    Y_train(f) = mean(correspondingRows);
end
X_train = C;


X_test = X_train;

[success,lscales,sigma_f,sigma_y] = GP_train_anisotropic(X_train,Y_train);
for f = 1:nfactorstoconsider
    fprintf('Factor %s:\t length-scale: %f\n',...
            var_names{f},lscales(f));
end
[opt_mu_s,opt_cov_s] = GPfun_posterior_anisotropic(X_test,X_train,Y_train,lscales,sigma_f,sigma_y);
opt_sigmas = sqrt(diag(opt_cov_s));


