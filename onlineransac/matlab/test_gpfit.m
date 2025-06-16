% Assume you have vectors x and y from your scatter data
x = x(:);  % Ensure column vectors
y = y(:);

% Fit GPR model using a squared exponential kernel
gprMdl = fitrgp(x, y, ...
    'KernelFunction', 'squaredexponential', ...
    'Sigma', std(y)*0.1, ...
    'Standardize', true);

% Predict using the trained model
xq = linspace(min(x), max(x), 500)';
[y_pred, y_sd] = predict(gprMdl, xq);

% Plot results
figure;
hold on;
scatter(x, y, 10, 'b', 'filled');         % Original data
plot(xq, y_pred, 'r-', 'LineWidth', 2);   % GPR mean prediction
fill([xq; flipud(xq)], ...
     [y_pred + 2*y_sd; flipud(y_pred - 2*y_sd)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Confidence band
xlabel('x'); ylabel('y');
title('Gaussian Process Regression Fit');
legend('Data', 'GPR Mean', '95% Confidence Interval');
hold off;
