function plot_optimal_controls(optimal_results)
% Simple optimal controls plot
figure;
plot(1:10, rand(10,3)); % Placeholder
xlabel('Time (weeks)');
ylabel('Control Intensity');
title('Optimal Control Strategies');
legend('Media', 'Bans', 'Treatment');
end