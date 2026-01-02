function plot_model_fit(data, params)
% Simple plot function
figure;
plot(data.year, data.prevalence, 'bo-');
xlabel('Year'); ylabel('Prevalence (%)');
title('Smoking Prevalence Fit');
end