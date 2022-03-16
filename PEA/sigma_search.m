function [mean_REE,coef] = sigma_search(sigma_z, rho_z , T, N_burnout, param, steady_state, quad)


%% parameter setup
alpha = param.alpha;
nss = steady_state.nss;
kss = steady_state.kss;
gamma = param.gamma;
mu = param.mu;
gx = param.gx;
dss = param.dss;
beta_s = param.beta_s;
phi = param.phi;
B = param.B;

% exogenous shock
sigmaz   = sigma_z;
rhoz = rho_z;
sigmat = param.sigmat;
rhot = param.rhot;



%% Upper bound and lower bound
N(1)=(nss)*0.5;  % the lower bound on n(t)
N(100)=(nss)*2;  % the upper bound on n(t)


%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1) = kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
u = zeros(T,1); % u(t)


coef = [log((1-nss)^(-mu)*(gx/beta_s)); 0.001; 0.001; 0.001]; % initial coefficients for labor parametrization

criter  	= 1e-6;            				        % Convergence criterion
update  	= .6;            				        % Updating rate 
%% Quadrature nodes, weight
n_nodes = quad.n_nodes;
epsi_nodes = quad.epsi_nodes;
weight_nodes = quad.weight_nodes;


%% Loop
iteration = 0; % initializing iteration count
dif = Inf;  % norm(zeta - coef)

while (dif > criter)
    up_bound  = kss * 1.5;     % Upper bound 150% of the steady state kss.
    low_bound = kss * 0.5;         % Lower bound 50% of the steady state kss.
for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) *kp(t) +  coef(3) * theta(t) + coef(4) * z(t)))^(-1/mu);
    
    n(t)=n(t)*(n(t)>=N(1))*(n(t)<=N(100))+N(1)*(n(t)<N(1))+N(100)*(n(t)>N(100)); %making n(t) not go over the bounds 
    % calculate u(t) using FOC
    u(t) = (alpha * theta(t) * kp(t)^(alpha-1) * n(t)^(1-alpha) / (z(t)* dss) )^(1/(phi-alpha));
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) * (1/((1-alpha)*theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(-alpha))))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (z(t)* dss/phi * u(t)^(phi)))*kp(t) + theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(1-alpha) - c(t));

    %calculate u(t) counterpart
    if kp(t+1) > up_bound  
          kp(t+1) = up_bound; 
   elseif kp(t+1) < low_bound
          kp(t+1) = low_bound;
    end

end



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

% Expectation part
e = (1-n_p).^(-mu) .* (theta.* (u.* k).^(alpha).* n.^(-alpha)) ./ (theta_p.*(u_p.*kprime).^(alpha).*n_p.^(-alpha)) .*...
    ((1- (z_p.*dss/phi .* u_p.^(phi))) + alpha .* theta_p .* u_p.^(alpha) .* kprime.^(alpha-1) .* n_p.^(1-alpha)) ;
e = e(1:end-1,:) *weight_nodes;

X = [ones(T-1,1), kp(1:end-2,:), theta(1:end-1,:), z(1:end-1,:)];
zeta = nlinfit(X,e,'object',coef);
dif = norm(coef - zeta);

if rem(iteration,100) == 0
    dif;
end
if dif > criter
    coef = update*zeta + (1-update)*coef;
%else
   % fprintf("Optimal Coefficient for parametrized function is %.3f \n", coef);
end
iteration = iteration+1;

end



%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1,1) = kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
u = zeros(T,1); % u(t)

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


kp = kp(N_burnout+1:end,:);
c = c(N_burnout+1:end,:);
z = z(N_burnout+1:end,:);
theta=theta(N_burnout+1:end,:);




%% Expectation part calculation
% z(t+1) 
z_p = repmat(z.^(rhoz),1,n_nodes) .* exp(sigmaz * epsi_nodes(:,2))';
% theta(t+1)
theta_p = repmat(theta.^(rhot),1,n_nodes) .* exp(sigmat * epsi_nodes(:,1))';

kprime = repmat(kp(2:end),1,n_nodes);
%n(t+1) - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* kprime +  coef(3) .* theta_p + coef(4).* z_p)).^(-1/mu));
%u(t+1) - FOC
u_p = (alpha .* theta_p .* kprime.^(alpha-1) .* n_p.^(1-alpha) ./ (z_p.* dss)).^(1/(phi-alpha));
%c(t+1) - FOC
c_p = ((B.*(1-n_p).^(-mu) ./ ((1-alpha)*theta_p.*(u_p.*kprime).^(alpha).*n_p.^(-alpha))).^(-1/gamma));





REE =  ((param.beta_s .* c_p.^(-param.gamma) .* (1- (param.dss/param.phi .* z_p .* u_p.^(param.phi))) + param.alpha .* theta_p .* u_p.^(param.alpha)...
    .* kprime.^(param.alpha-1) .* n_p.^(1-param.alpha)) ./ (param.gx .* c.^(-param.gamma)))*weight_nodes  - 1;

mean_REE = mean(log10(abs(REE)));




