function [c, kp, n, u, y, w,i]= simulation(sigma_z, rho_z, coef , T, N_burnout, param, steady_state)
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
y = zeros(T,1); % y(t)
w = zeros(T,1); % w(t)
i = zeros(T,1);

for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) *kp(t) +  coef(3) * theta(t) + coef(4) * z(t)))^(-1/mu);;
    % calculate u(t) using FOC
    u(t) = (alpha * theta(t) * kp(t)^(alpha-1) * n(t)^(1-alpha) / (z(t)* dss) )^(1/(phi-alpha));
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) * (1/((1-alpha)*theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(-alpha))))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (z(t)* dss/phi * u(t)^(phi)))*kp(t) + theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(1-alpha) - c(t));
    %calculate output
    y(t)= theta(t) * kp(t)^(alpha)*n(t)^(1-alpha)*u(t)^(alpha);
    %calculate marginal product of labor
    w(t)= (1-alpha) * theta(t) * (u(t)* kp(t))^(alpha) * n(t)^(-alpha);
    %calculate investment
    i(t) = y(t) - c(t);
    

end


i = i(N_burnout+1:end,:);
n = n(N_burnout+1:end,:);
u = u(N_burnout+1:end,:);
kp = kp(N_burnout+1:end,:);
c = c(N_burnout+1:end,:);
y = y(N_burnout+1:end,:);
w = w(N_burnout+1:end,:);

