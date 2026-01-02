function R_s = calculate_R0(params)
% Calculate basic reproduction number
a = params.alpha1 + params.alpha2 + params.alpha3 + params.mu + params.d1;
b = params.mu + params.d2 + params.gamma;
R_s = (params.alpha2 * params.beta) / (a * b);
end