function loglike = likelihood(data, params)
% Likelihood function for Bayesian estimation
% Compares model predictions to observed prevalence data

% Simulate model
[t, Y] = simulate_model(params, 2000, 2022);
N = sum(Y(:,1:6), 2);
model_prevalence = (Y(:,2) + Y(:,3)) ./ N * 100;

% Interpolate to match data years
model_years = 2000:2022;
data_years = data.year;

% Find common years
[common_years, idx_model, idx_data] = intersect(model_years, data_years);

if isempty(common_years)
    loglike = -inf;
    return;
end

% Calculate likelihood (Gaussian errors)
model_vals = model_prevalence(idx_model);
data_vals = data.prevalence(idx_data);
errors = model_vals - data_vals;

% Assume sigma = 0.5% measurement error
sigma = 0.5;
loglike = sum(-0.5 * (errors/sigma).^2 - log(sigma*sqrt(2*pi)));

% Penalize if model prevalence is negative
if any(model_vals < 0)
    loglike = loglike - 1000;
end
end