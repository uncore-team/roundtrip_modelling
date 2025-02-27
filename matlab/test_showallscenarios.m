% show basically all scenarios
%
% We can separate scenarios #7-#17 by using their std: they have larger std
% than the rest, except #11 and #13. #11 has larger skewness, but #13 has
% no differetiation; anyway, both have similar behavior as the rest not
% included in #7-#17, thus we conclude that a large std produces quite
% different effect in the difference between the expectation of the current 
% model and the expectation of the future model.

close all;
clear;
clc;

figscssep = 0;  % 1 for separate figures for scenarios; 0 for one fig

expcat = ExperimentCatalog(1);
numexps = length(expcat);
indstartsc = 1;
indendsc = 200;

hs1 = nan(1,numexps);
hs2 = nan(1,numexps);
hs3 = nan(1,numexps);
labs = cell(1,numexps);
colors = 'kbrmcg';
styles = '.*os^v';
styles2 = {'-','--',':','-.','none','none'};
lensts = length(styles);
if lensts ~= length(colors)
    error('Unmatched colors and styles');
end


netfact = nan(1,numexps);
densfact = nan(1,numexps);
distfact = nan(1,numexps);
clfact = nan(1,numexps);
sefact = nan(1,numexps);
for f = 1:numexps
    [~,levsnum] = ExperimentFactors(expcat{f}.class,expcat{f}.index,...
                                    'bernoulli',[]); % algorithm is irrelevant for this loop
    distfact(f) = levsnum(1); % geographical distance
    sefact(f) = levsnum(2);   % server sw
    clfact(f) = levsnum(3);   % client sw
    netfact(f) = levsnum(4);  % network
    densfact(f) = levsnum(5); % density
end
data = [distfact.' sefact.' clfact.' netfact.' densfact.'];
numVars = size(data,2);
varNames = {'dist', 'server', 'client', 'network', 'density'};


% Normalize the data to range [0, 1] for better comparison
minData = min(data, [], 1);
maxData = max(data, [], 1);
normData = (data - minData) ./ (maxData - minData);
% Angles for the radar chart
theta = linspace(0, 2*pi, numVars + 1);
% Create figure
figure(101);
hold on;
% Plot each data point
for i = 1:size(normData, 1)
    r = normData(i, :);
    r = [r, r(1)]; % Close the loop
    x = r .* cos(theta);
    y = r .* sin(theta);
    plot(x, y, '-o', 'LineWidth', 2);
end
% Draw radial lines and labels
for i = 1:numVars
    angle = theta(i);
    x = [0, cos(angle)];
    y = [0, sin(angle)];
    plot(x, y, 'k-'); % Radial lines
    text(cos(angle), sin(angle), sprintf('%s', varNames{i}), 'HorizontalAlignment', 'center');
end
% Set axis limits
axis equal;
axis([-1 1 -1 1]);
% Add title
title('Radar Chart of 5-Dimensional Data');
% Release hold
hold off;

% data_mean = mean(data);
% data_centered = data - data_mean;
% cov_matrix = cov(data_centered);
% % Perform eigenvalue decomposition on the covariance matrix
% [eigenvectors, eigenvalues_matrix] = eig(cov_matrix);
% % Extract the eigenvalues from the diagonal matrix
% eigenvalues = diag(eigenvalues_matrix);
% % Sort the eigenvalues and corresponding eigenvectors in descending order
% [eigenvalues_sorted, sort_indices] = sort(eigenvalues, 'descend');
% eigenvectors_sorted = eigenvectors(:, sort_indices);
% % Calculate the principal components
% principal_components = data_centered * eigenvectors_sorted;
% % Plot the first two principal components (for visualization purposes)
% figure;
% scatter(principal_components(:, 1), principal_components(:, 2), 50, 'filled');
% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% title('PCA Result');
% grid on;
% % Display the percentage of variance explained by each principal component
% variance_explained = eigenvalues_sorted / sum(eigenvalues_sorted) * 100;
% disp('Variance explained by each principal component:');
% disp(variance_explained);
% % Optionally, you can choose the number of principal components to keep
% % For example, keep the first 2 principal components
% num_components_to_keep = 2;
% data_reduced = principal_components(:, 1:num_components_to_keep);



% ---------- Figures

figure(100);
[h, ax] = plotmatrix(data);
for i = 1:numVars
    for j = 1:numVars
        if j == 1
            ylabel(ax(i,j), varNames{i}, 'FontWeight', 'bold');
        end
        if i == numVars
            xlabel(ax(i,j), varNames{j}, 'FontWeight', 'bold');
        end
        grid(ax(i,j), 'on');
    end
end

figure(1);
hold on;
grid;
figure(2);
hold on;
grid;
figure(3);
subplot(1,3,1);
hold on;
grid;
subplot(1,3,2);
hold on;
grid;
subplot(1,3,3);
hold on;
grid;
meansscs = zeros(1,numexps);
for f = 1:numexps
    [~,~,data] = ExperimentGet(expcat{f}.class,expcat{f}.index,1,Inf,0,NaN,0);
    data = data(indstartsc:indendsc);
    ts = cumsum(data);
    meansscs(f) = mean(data);
    indcols = mod(f-1,lensts)+1;
    indstyls = floor((f-1)/lensts)+1;
    commcol = sprintf('%c-%c',colors(indcols),styles(indstyls));
    commcol3 = sprintf('%c%c',colors(indcols),styles(indstyls));
    
    if figscssep
        figure;
    else
        figure(1);
    end
    hs1(f) = plot(ts,data,commcol);

    figure(2);
    validdata = data((~isnan(data)) & (~isinf(data)));
    mu = mean(validdata);
    si = std(validdata);
    sk = skewness(validdata);
    if ~isempty(validdata)
        hs2(f) = histogram(validdata,'Normalization','pdf','FaceColor',colors(indcols),'LineStyle',styles2(indstyls));
    else
        hs2(f) = plot(0,0,commcol);
    end

    figure(3);
    subplot(1,3,1);
    hs3(f) = plot(mu,f,commcol3);
    subplot(1,3,2);
    plot(si,f,commcol3);
    subplot(1,3,3);
    plot(sk,f,commcol3);

    labs{f} = sprintf('#%d',f);
end
figure(1);
xlabel('Time (ms)');
ylabel('Delay (ms)');
legend(hs1,labs);
figure(2);
xlabel('Delay (ms)');
ylabel('Density');
legend(hs2,labs);
figure(3);
subplot(1,3,1);
legend(hs3,labs);
ylabel('exp (#)');
xlabel('Avg (ms)');
subplot(1,3,2);
ylabel('exp (#)');
xlabel('Std (ms)');
subplot(1,3,3);
ylabel('exp (#)');
xlabel('Skewness');
