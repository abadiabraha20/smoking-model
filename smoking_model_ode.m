function dYdt = smoking_model_ode(t, Y, params)
% ODE system for 7-compartment smoking model
% Y = [P, L, H, Qt, Qp, Sm, M]

% Unpack parameters
beta = params.beta;
alpha1 = params.alpha1;
alpha2 = params.alpha2;
alpha3 = params.alpha3;
gamma = params.gamma;
epsilon0 = params.epsilon0;
mu = params.mu;
d1 = params.d1;
d2 = params.d2;
d3 = params.d3;
sigma = params.sigma;
eta = params.eta;
varphi = params.varphi;
phi = params.phi;

% Unpack states
P = Y(1); L = Y(2); H = Y(3); Qt = Y(4); Qp = Y(5); Sm = Y(6); M = Y(7);

% Total human population
N = P + L + H + Qt + Qp + Sm;

% System of differential equations
dPdt = mu*N - (beta*P*H)/N - epsilon0*P*M - mu*P;
dLdt = (beta*P*H)/N - (alpha1 + alpha2 + alpha3 + mu + d1)*L;
dHdt = alpha2*L - (mu + d2 + gamma)*H;
dQtdt = alpha3*L + gamma*sigma*H - (eta + mu + d3)*Qt;
dQpdt = alpha1*L + gamma*(1-sigma)*H + eta*Qt - mu*Qp;
dSmdt = epsilon0*P*M - (varphi*L + mu)*Sm;
dMdt = varphi*L - phi*M;

% Return derivatives
dYdt = [dPdt; dLdt; dHdt; dQtdt; dQpdt; dSmdt; dMdt];
end