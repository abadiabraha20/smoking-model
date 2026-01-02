function prcc_results = calculate_prcc(sensitivity_results)
% Calculate Partial Rank Correlation Coefficients (PRCC)

fprintf('Calculating PRCC...\n');

param_names = sensitivity_results.param_names;
n_params = length(param_names);
n_outputs = 4; % R_s, endemic, peak, time

% Output matrices
outputs = [sensitivity_results.R_s, ...
           sensitivity_results.endemic_prevalence, ...
           sensitivity_results.peak_prevalence, ...
           sensitivity_results.time_to_peak];
output_names = {'R_s', 'Endemic_Prevalence', 'Peak_Prevalence', 'Time_to_Peak'};

% Initialize PRCC storage
prcc_results.PRCC = zeros(n_params, n_outputs);
prcc_results.p_values = zeros(n_params, n_outputs);

% Calculate PRCC for each output
for out_idx = 1:n_outputs
    Y = outputs(:, out_idx);
    
    for param_idx = 1:n_params
        X = sensitivity_results.parameters(:, param_idx);
        
        % Rank transform
        X_rank = tiedrank(X);
        Y_rank = tiedrank(Y);
        
        % Calculate partial correlation
        other_params = sensitivity_results.parameters;
        other_params(:, param_idx) = []; % Remove current parameter
        
        if ~isempty(other_params)
            % Rank transform other parameters
            other_ranks = tiedrank(other_params);
            
            % Residuals after removing effect of other parameters
            X_resid = X_rank - other_ranks * ...
                     (other_ranks \ X_rank);
            Y_resid = Y_rank - other_ranks * ...
                     (other_ranks \ Y_rank);
            
            % Correlation of residuals
            prcc = corr(X_resid, Y_resid);
        else
            prcc = corr(X_rank, Y_rank);
        end
        
        % Significance test
        n = length(Y);
        t_stat = prcc * sqrt((n-2) / (1-prcc^2));
        p_value = 2 * (1 - tcdf(abs(t_stat), n-2));
        
        prcc_results.PRCC(param_idx, out_idx) = prcc;
        prcc_results.p_values(param_idx, out_idx) = p_value;
    end
end

% Store metadata
prcc_results.param_names = param_names;
prcc_results.output_names = output_names;

% Display top parameters
fprintf('\n=== PRCC Results (|PRCC| > 0.3) ===\n');
for out_idx = 1:n_outputs
    fprintf('\n%s:\n', output_names{out_idx});
    [sorted_prcc, sort_idx] = sort(abs(prcc_results.PRCC(:, out_idx)), 'descend');
    
    for i = 1:min(5, n_params)
        if sorted_prcc(i) > 0.3
            param_idx = sort_idx(i);
            fprintf('  %-12s: PRCC = %6.3f (p = %.4f)\n', ...
                    param_names{param_idx}, ...
                    prcc_results.PRCC(param_idx, out_idx), ...
                    prcc_results.p_values(param_idx, out_idx));
        end
    end
end
end