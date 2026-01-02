function save_results(params, sensitivity_results, optimal_results)
% Save results
if ~exist('Results', 'dir')
    mkdir('Results');
end
save('Results/analysis_results.mat', 'params', 'sensitivity_results', 'optimal_results');
fprintf('Results saved to Results/ folder\n');
end