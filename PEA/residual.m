function [mean_REE, max_REE, mean_RBC, max_RBC, mean_RMUL, max_RMUL] = residual(sigma_z, rho_z, coef , T, N_burnout, param, steady_state, quad)
%% parameter setup
alpha = param.alpha;
kss = steady_state.kss;
gamma = param.gamma;
mu = param.mu;
gx = param.gx;

sigmaz   = sigma_z;
rhoz = rho_z;
sigmat = param.sigmat;
rhot = param.rhot;

dss = param.dss;
beta_s = param.beta_s;
phi = param.phi;
B = param.B;

%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1,1) = kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
u = zeros(T,1); % u(t)


%% Quadrature
n_nodes = quad.n_nodes;
epsi_nodes = quad.epsi_nodes;
weight_nodes = quad.weight_nodes;

for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) *kp(t) +  coef(3) * theta(t) + coef(4) * z(t)))^(-1/mu);
    % calculate u(t) using FOC
    u(t) = (alpha * theta(t) * kp(t)^(alpha-1) * n(t)^(1-alpha) / (z(t)* dss) )^(1/(phi-alpha));
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) * (1/((1-alpha)*theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(-alpha))))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (z(t)* dss/phi * u(t)^(phi)))*kp(t) + theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(1-alpha) - c(t));

end



n = n(N_burnout+1:end,:);
u = u(N_burnout+1:end,:);
kp = kp(N_burnout+1:end,:);
c = c(N_burnout+1:end,:);
z = z(N_burnout+1:end,:);
theta=theta(N_burnout+1:end,:);




%% Expectation part calculation
% z(t+1) 
z_p = repmat(z.^(rhoz),1,n_nodes) .* exp(sigmaz * epsi_nodes(:,2))';
% theta(t+1)
theta_p = repmat(theta.^(rhot),1,n_nodes) .* exp(sigmat * epsi_nodes(:,1))';
%k(t) and k(t+1)
k = repmat(kp(1:end-1),1,n_nodes);
kprime = repmat(kp(2:end),1,n_nodes);
%n(t+1) - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* kprime +  coef(3) .* theta_p + coef(4).* z_p)).^(-1/mu));
%u(t+1) - FOC
u_p = (alpha .* theta_p .* kprime.^(alpha-1) .* n_p.^(1-alpha) ./ (z_p.* dss)).^(1/(phi-alpha));
%c(t+1) - FOC
c_p = ((B.*(1-n_p).^(-mu) ./ ((1-alpha)*theta_p.*(u_p.*kprime).^(alpha).*n_p.^(-alpha))).^(-1/gamma));


%% compute residuals 
uniform_weight = ones(n_nodes,1)./n_nodes;
RBC = ((1- (z.* param.dss/param.phi .* u.^(param.phi))) .* k + theta .* (u .* k).^(param.alpha) .* n.^(1-param.alpha)) ./ (c + param.gx .* kprime(:,1)) * uniform_weight  -  1;

REE =  ((param.beta_s .* c_p.^(-param.gamma) .* (1- (param.dss/param.phi .* z_p .* u_p.^(param.phi))) + param.alpha .* theta_p .* u_p.^(param.alpha)...
    .* kprime.^(param.alpha-1) .* n_p.^(1-param.alpha)) ./ (param.gx .* c.^(-param.gamma)))*weight_nodes  - 1;

RMUL = ((c.^(-param.gamma) .* (1-param.alpha).*theta.*(u.*k).^(param.alpha).*n.^(-param.alpha))./ (param.B.*(1-n).^(-param.mu))) * uniform_weight  -1;


max_RBC = max(log10(abs(RBC)));
mean_RBC = mean(log10(abs(RBC)));

max_REE = max(log10(abs(REE)));
mean_REE = mean(log10(abs(REE)));

max_RMUL = max(log10(abs(RMUL)));
mean_RMUL = mean(log10(abs(RMUL)));



