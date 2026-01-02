%% MAIN SCRIPT: Smoking Dynamics Model Analysis
% Reproduces all results from the manuscript

clear all; close all; clc;
fprintf('=== Smoking Dynamics Model Analysis ===\n');
fprintf('Starting at %s\n', datestr(now));

%% Add paths
addpath(genpath('1_Model_Simulation'));
addpath(genpath('2_Bayesian_Estimation'));
addpath(genpath('3_Sensitivity_Analysis'));
addpath(genpath('4_Optimal_Control'));
addpath(genpath('5_Figure_Generation'));
addpath(genpath('6_Data_Management'));
addpath(genpath('7_Utilities'));

%% Step 1: Load data
fprintf('\n1. Loading Ethiopian smoking data...\n');
data = load_data();
years = data.year;
prevalence = data.prevalence;

%% Step 2: Set default parameters
fprintf('2. Setting model parameters...\n');
params = default_parameters();

%% Step 3: Run baseline simulation
fprintf('3. Running baseline simulation...\n');
[t, Y] = simulate_model(params, 0, 200); % 200 weeks
N = sum(Y(:,1:6), 2);
baseline_prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
fprintf('   Endemic prevalence: %.1f%%\n', baseline_prevalence(end));

%% Step 4: Calculate reproduction number
fprintf('4. Calculating reproduction number...\n');
R_s = calculate_R0(params);
fprintf('   R_s = %.3f\n', R_s);

%% Step 5: Find equilibria
fprintf('5. Finding equilibria...\n');
[DFE, EE] = find_equilibria(params);
fprintf('   Smoking-free equilibrium found\n');
fprintf('   Endemic equilibrium prevalence: %.1f%%\n', EE.prevalence*100);

%% Step 6: Sensitivity analysis
fprintf('6. Running sensitivity analysis...\n');
sensitivity_results = lhs_sampling(params, 1000); % 1000 samples for speed
prcc_results = calculate_prcc(sensitivity_results);

%% Step 7: Optimal control
fprintf('7. Solving optimal control problem...\n');
optimal_results = solve_optimal_control(params, 52); % 52 weeks
fprintf('   Prevalence reduction: %.1f%%\n', optimal_results.prevalence_reduction);

%% Step 8: Generate key figures
fprintf('8. Generating manuscript figures...\n');
plot_model_fit(data, params);
plot_sensitivity_results(prcc_results);
plot_optimal_controls(optimal_results);

%% Step 9: Save results
fprintf('9. Saving results...\n');
save_results(params, sensitivity_results, optimal_results);

fprintf('\n=== Analysis complete ===\n');
fprintf('Results saved in Results/ folder\n');
fprintf('Figures saved in Figures/ folder\n');