function sobol_results = sobol_indices(params, n_samples)
% Calculate Sobol sensitivity indices

fprintf('Calculating Sobol indices (N=%d)...\n', n_samples);

% Parameter ranges
param_ranges = {
    'beta',     [0.40, 0.70];
    'alpha2',   [0.65, 1.00];
    'gamma',    [0.20, 0.40];
    'epsilon0', [0.08, 0.16];
};
n_params = length(param_ranges);

% Generate samples for Sobol analysis
A = lhsdesign(n_samples/2, n_params);
B = lhsdesign(n_samples/2, n_params);

% Initialize
Y_A = zeros(n_samples/2, 1);
Y_B = zeros(n_samples/2, 1);
Y_C = zeros(n_params, n_samples/2);

% Evaluate model for sample sets
fprintf('  Evaluating model...\n');
for i = 1:n_samples/2
    % Sample set A
    params_A = params;
    for j = 1:n_params
        param_name = param_ranges{j,1};
        min_val = param_ranges{j,2}(1);
        max_val = param_ranges{j,2}(2);
        params_A.(param_name) = min_val + A(i,j) * (max_val - min_val);
    end
    [~, Y] = simulate_model(params_A, 0, 200);
    N = sum(Y(:,1:6), 2);
    prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
    Y_A(i) = mean(prevalence(end-50:end));
    
    % Sample set B
    params_B = params;
    for j = 1:n_params
        param_name = param_ranges{j,1};
        min_val = param_ranges{j,2}(1);
        max_val = param_ranges{j,2}(2);
        params_B.(param_name) = min_val + B(i,j) * (max_val - min_val);
    end
    [~, Y] = simulate_model(params_B, 0, 200);
    N = sum(Y(:,1:6), 2);
    prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
    Y_B(i) = mean(prevalence(end-50:end));
    
    % Sample sets C_i (only parameter i from B, others from A)
    for j = 1:n_params
        params_C = params;
        for k = 1:n_params
            param_name = param_ranges{k,1};
            min_val = param_ranges{k,2}(1);
            max_val = param_ranges{k,2}(2);
            if k == j
                value = min_val + B(i,k) * (max_val - min_val);
            else
                value = min_val + A(i,k) * (max_val - min_val);
            end
            params_C.(param_name) = value;
        end
        [~, Y] = simulate_model(params_C, 0, 200);
        N = sum(Y(:,1:6), 2);
        prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;
        Y_C(j,i) = mean(prevalence(end-50:end));
    end
    
    if mod(i, 100) == 0
        fprintf('    Sample %d/%d\n', i, n_samples/2);
    end
end

% Calculate Sobol indices
f0 = mean(Y_A);
D = mean(Y_A.^2) - f0^2;

first_order = zeros(n_params, 1);
total_order = zeros(n_params, 1);

for i = 1:n_params
    % First-order index (main effect)
    first_order(i) = (mean(Y_A .* Y_C(i,:)') - f0^2) / D;
    
    % Total-effect index
    total_order(i) = 1 - (mean(Y_B .* Y_C(i,:)') - f0^2) / D;
end

% Store results
sobol_results.first_order = first_order;
sobol_results.total_order = total_order;
sobol_results.param_names = param_ranges(:,1);
sobol_results.interaction = total_order - first_order;

% Display
fprintf('\n=== Sobol Indices ===\n');
fprintf('%-12s %-12s %-12s %-12s\n', 'Parameter', 'First-Order', 'Total-Order', 'Interaction');
fprintf('%s\n', repmat('-', 1, 48));
for i = 1:n_params
    fprintf('%-12s %-12.4f %-12.4f %-12.4f\n', ...
            param_ranges{i,1}, ...
            first_order(i), ...
            total_order(i), ...
            sobol_results.interaction(i));
end
end