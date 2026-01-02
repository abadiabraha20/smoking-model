function plot_posteriors(results)
% Plot posterior distributions from MCMC

param_names = results.param_names;
n_params = length(param_names);

figure('Position', [100, 100, 1200, 800]);

for i = 1:n_params
    subplot(2, ceil(n_params/2), i);
    
    chain = results.chains(:, i);
    param_name = param_names{i};
    
    % Histogram
    histogram(chain, 30, 'Normalization', 'pdf', ...
              'FaceColor', [0.2, 0.4, 0.6], 'EdgeColor', 'none');
    hold on;
    
    % Add mean and CI
    mean_val = mean(chain);
    CI = quantile(chain, [0.025, 0.975]);
    
    plot([mean_val, mean_val], ylim, 'r-', 'LineWidth', 2);
    plot([CI(1), CI(1)], ylim, 'r--', 'LineWidth', 1.5);
    plot([CI(2), CI(2)], ylim, 'r--', 'LineWidth', 1.5);
    
    % Labels
    xlabel(param_name, 'FontSize', 11);
    ylabel('Density', 'FontSize', 11);
    title(sprintf('%s: %.3f [%.3f, %.3f]', ...
          param_name, mean_val, CI(1), CI(2)), 'FontSize', 10);
    grid on;
    box on;
end

sgtitle('Posterior Distributions of Estimated Parameters', 'FontSize', 14);

% Save figure
if ~exist('Figures', 'dir')
    mkdir('Figures');
end
saveas(gcf, 'Figures/posterior_distributions.png');
saveas(gcf, 'Figures/posterior_distributions.fig');
end