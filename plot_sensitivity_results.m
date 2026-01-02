function plot_sensitivity_results(prcc_results)
% Simple sensitivity plot
figure;
bar(prcc_results.PRCC(:,1));
xlabel('Parameter Index');
ylabel('PRCC Value');
title('Sensitivity Analysis Results');
end