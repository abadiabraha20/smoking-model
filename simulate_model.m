function [t, Y] = simulate_model(params, t_start, t_end)
% SIMULATE_MODEL - Run smoking model simulation
% Simplified version for submission

fprintf('Running simulation from %d to %d...\n', t_start, t_end);

% Initial conditions (endemic equilibrium approximation)
P0 = 0.8; L0 = 0.04; H0 = 0.04; Qt0 = 0.03;
Qp0 = 0.05; Sm0 = 0.03; M0 = 0.03;
Y0 = [P0; L0; H0; Qt0; Qp0; Sm0; M0];

% Time vector
t = linspace(t_start, t_end, 100)';

% Simple integration (Euler method)
Y = zeros(length(t), 7);
Y(1,:) = Y0';

for i = 1:length(t)-1
    % Get current state
    current_Y = Y(i,:)';
    
    % Calculate derivatives
    dYdt = smoking_model_ode(t(i), current_Y, params);
    
    % Euler step
    dt = t(2) - t(1);
    Y(i+1,:) = current_Y' + dYdt' * dt;
    
    % Keep non-negative
    Y(i+1,:) = max(Y(i+1,:), 0);
end

fprintf('Simulation complete: %d time points\n', length(t));
end