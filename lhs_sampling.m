function results = lhs_sampling(params, n_samples)
% Latin Hypercube Sampling for global sensitivity analysis

fprintf('Running LHS with %d samples...\n', n_samples);

% Parameter ranges for sensitivity analysis
param_ranges = {
    'beta',     [0.40, 0.70];
    'alpha1',   [0.06, 0.14];
    'alpha2',   [0.65, 1.00];
    'alpha3',   [0.20, 0.40];
    'gamma',    [0.20, 0.40];
    'epsilon0', [0.08, 0.16];
    'mu',       [0.006, 0.009];
    'd1',       [0.004, 0.010];
    'd2',       [0.002, 0.006];
    'sigma',    [0.40, 0.60];
    'varphi',   [0.0004, 0.0010];
    'phi',      [0.05, 0.08];
};

n_params = length(param_ranges);

% Generate LHS design
lhs_design = lhsdesign(n_samples, n_params);

% Initialize storage
results.parameters = zeros(n_samples, n_params);
results.R_s = zeros(n_samples, 1);
results.endemic_prevalence = zeros(n_samples, 1);
results.peak_prevalence = zeros(n_samples, 1);
results.time_to_peak = zeros(n_samples, 1);
results.param_names = param_ranges(:,1);

% Sample and evaluate model
for i = 1:n_samples
    % Create parameter set
    sample_params = params;
    for j = 1:n_params
        param_name = param_ranges{j,1};
        min_val = param_ranges{j,2}(1);
        max_val = param_ranges{j,2}(2);
        
        value = min_val + lhs_design(i,j) * (max_val - min_val);
        sample_params.(param_name) = value;
        results.parameters(i,j) = value;
    end
    
    % Calculate R_s
    results.R_s(i) = calculate_R0(sample_params);
    
    % Simulate to find endemic prevalence
    [t, Y] = simulate_model(sample_params, 0, 500);
    N = sum(Y(:,1:6), 2);
    prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
    
    results.endemic_prevalence(i) = mean(prevalence(end-100:end));
    results.peak_prevalence(i) = max(prevalence);
    
    % Find time to peak
    [~, peak_idx] = max(prevalence);
    results.time_to_peak(i) = t(peak_idx);
    
    % Progress
    if mod(i, 200) == 0
        fprintf('  Completed %d/%d samples\n', i, n_samples);
    end
end

fprintf('LHS sampling complete.\n');
end