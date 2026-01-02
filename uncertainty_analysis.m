function uncertainty_results = uncertainty_analysis(params, n_samples)
% Uncertainty propagation through model

fprintf('Running uncertainty propagation (%d samples)...\n', n_samples);

% Generate parameter samples
param_ranges = {
    'beta',     [0.447, 0.617];  % 95% CI from Bayesian
    'alpha2',   [0.698, 0.932];
    'gamma',    [0.257, 0.353];
    'epsilon0', [0.099, 0.149];
};
n_params = length(param_ranges);

samples = zeros(n_samples, n_params);
for i = 1:n_params
    min_val = param_ranges{i,2}(1);
    max_val = param_ranges{i,2}(2);
    samples(:,i) = min_val + (max_val - min_val) * rand(n_samples, 1);
end

% Run model for all samples
R_s_samples = zeros(n_samples, 1);
prevalence_samples = zeros(n_samples, 1);

for s = 1:n_samples
    % Create parameter set
    sample_params = params;
    for i = 1:n_params
        param_name = param_ranges{i,1};
        sample_params.(param_name) = samples(s,i);
    end
    
    % Calculate R_s
    R_s_samples(s) = calculate_R0(sample_params);
    
    % Simulate to endemic equilibrium
    [~, Y] = simulate_model(sample_params, 0, 300);
    N = sum(Y(:,1:6), 2);
    prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
    prevalence_samples(s) = mean(prevalence(end-50:end));
    
    if mod(s, 200) == 0
        fprintf('  Sample %d/%d\n', s, n_samples);
    end
end

% Calculate statistics
uncertainty_results.R_s.mean = mean(R_s_samples);
uncertainty_results.R_s.median = median(R_s_samples);
uncertainty_results.R_s.CI = quantile(R_s_samples, [0.025, 0.975]);
uncertainty_results.R_s.prob_above_1 = mean(R_s_samples > 1);

uncertainty_results.prevalence.mean = mean(prevalence_samples);
uncertainty_results.prevalence.median = median(prevalence_samples);
uncertainty_results.prevalence.CI = quantile(prevalence_samples, [0.025, 0.975]);

% Display results
fprintf('\n=== Uncertainty Analysis Results ===\n');
fprintf('R_s: %.3f [%.3f, %.3f]\n', ...
        uncertainty_results.R_s.median, ...
        uncertainty_results.R_s.CI(1), ...
        uncertainty_results.R_s.CI(2));
fprintf('P(R_s > 1) = %.3f\n', uncertainty_results.R_s.prob_above_1);
fprintf('\nEndemic prevalence: %.2f%% [%.2f%%, %.2f%%]\n', ...
        uncertainty_results.prevalence.median, ...
        uncertainty_results.prevalence.CI(1), ...
        uncertainty_results.prevalence.CI(2));

% Plot distributions
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
histogram(R_s_samples, 30, 'Normalization', 'pdf', ...
          'FaceColor', [0.2, 0.4, 0.6]);
hold on;
plot([1, 1], ylim, 'r--', 'LineWidth', 2);
xlabel('R_s', 'FontSize', 12);
ylabel('Density', 'FontSize', 12);
title('Uncertainty in Reproduction Number', 'FontSize', 14);
grid on;
box on;

subplot(1,2,2);
histogram(prevalence_samples, 30, 'Normalization', 'pdf', ...
          'FaceColor', [0.6, 0.2, 0.4]);
xlabel('Endemic Prevalence (%)', 'FontSize', 12);
ylabel('Density', 'FontSize', 12);
title('Uncertainty in Endemic Prevalence', 'FontSize', 14);
grid on;
box on;

% Save figure
if ~exist('Figures', 'dir')
    mkdir('Figures');
end
saveas(gcf, 'Figures/uncertainty_analysis.png');
saveas(gcf, 'Figures/uncertainty_analysis.fig');
end