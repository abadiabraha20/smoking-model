function [DFE, EE] = find_equilibria(params)
% Find disease-free and endemic equilibria

% Disease-free equilibrium (trivial)
DFE.P = 1; DFE.L = 0; DFE.H = 0; DFE.Qt = 0;
DFE.Qp = 0; DFE.Sm = 0; DFE.M = 0;
DFE.prevalence = 0;

% Endemic equilibrium (approximate analytical solution)
a = params.alpha1 + params.alpha2 + params.alpha3 + params.mu + params.d1;
b = params.mu + params.d2 + params.gamma;
R_s = (params.alpha2 * params.beta) / (a * b);

if R_s > 1
    % Calculate endemic equilibrium components
    H_star = (a*b*(R_s - 1)) / (params.beta*(params.alpha2 + b));
    L_star = (b/params.alpha2) * H_star;
    P_star = (a*L_star) / (params.beta*H_star);
    
    % Other compartments
    Qt_star = (params.alpha3*L_star + params.gamma*params.sigma*H_star) / ...
              (params.eta + params.mu + params.d3);
    Qp_star = (params.alpha1*L_star + params.gamma*(1-params.sigma)*H_star + ...
               params.eta*Qt_star) / params.mu;
    M_star = (params.varphi/params.phi) * L_star;
    
    % Media-aware compartment
    Sm_star = (params.epsilon0*P_star*M_star) / (params.varphi*L_star + params.mu);
    
    % Normalize
    total = P_star + L_star + H_star + Qt_star + Qp_star + Sm_star;
    EE.P = P_star/total;
    EE.L = L_star/total;
    EE.H = H_star/total;
    EE.Qt = Qt_star/total;
    EE.Qp = Qp_star/total;
    EE.Sm = Sm_star/total;
    EE.M = M_star;
    EE.prevalence = (L_star + H_star)/total;
else
    EE = DFE; % No endemic equilibrium if R_s <= 1
end
end