function metrics = model_validation(data, params)
% Calculate model validation metrics

% Simulate model
[t, Y] = simulate_model(params, 2000, 2022);
N = sum(Y(:,1:6), 2);
model_prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;

% Match years
model_years = 2000:2022;
[common_years, idx_model, idx_data] = intersect(model_years, data.year);

if isempty(common_years)
    metrics = struct();
    return;
end

model_vals = model_prevalence(idx_model);
data_vals = data.prevalence(idx_data);

% Calculate metrics
metrics.RMSE = sqrt(mean((model_vals - data_vals).^2));
metrics.MAE = mean(abs(model_vals - data_vals));
metrics.R2 = 1 - sum((model_vals - data_vals).^2) / ...
                sum((data_vals - mean(data_vals)).^2);

% Display
fprintf('\n=== Model Validation Metrics ===\n');
fprintf('RMSE: %.3f%%\n', metrics.RMSE);
fprintf('MAE:  %.3f%%\n', metrics.MAE);
fprintf('R²:   %.3f\n', metrics.R2);

% Plot fit
figure('Position', [100, 100, 800, 600]);
plot(data.year, data.prevalence, 'bo', 'MarkerSize', 8, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Observed');
hold on;
plot(model_years, model_prevalence, 'r-', 'LineWidth', 2, ...
     'DisplayName', 'Model fit');
xlabel('Year', 'FontSize', 12);
ylabel('Smoking Prevalence (%)', 'FontSize', 12);
title(sprintf('Model Validation (R² = %.3f)', metrics.R2), 'FontSize', 14);
legend('Location', 'best');
grid on;
box on;

% Save figure
if ~exist('Figures', 'dir')
    mkdir('Figures');
end
saveas(gcf, 'Figures/model_validation.png');
saveas(gcf, 'Figures/model_validation.fig');
end