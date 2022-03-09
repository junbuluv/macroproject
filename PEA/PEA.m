function mean_REE = PEA(input, T, param, steady_state, quad)
%% parameter setup
alpha = param.alpha;
nss = steady_state.nss;
kss = steady_state.kss;
gamma = param.gamma;
mu = param.mu;
gx = param.gx;


sigmaz   = input(1);
rhoz = input(2);
sigmat = param.sigmat;
rhot = param.rhot;

dss = param.dss;
beta_s = param.beta_s;
phi1 = param.phi1;
phi2 = param.phi2;
B = param.B;


%% Upper bound and lower bound
N(1)=(nss)*0.5;  % the lower bound on n(t)
N(100)=(nss)*2;  % the upper bound on n(t)


%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1) = 0.8*kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
e = zeros(T-1,1); % Parameterized expectation 
z_p = zeros(T-1,1); % z_p(t)
c_p = zeros(T-1,1); % c_p(t)
u = zeros(T,1); % u(t)
u_p = zeros(T-1,1);% u_p(t)
theta_p = zeros(T-1,1); % theta_p(t) 
n_p = zeros(T-1,1); %n_p(t)
kprime = zeros(T,1);


coef = [log((1-nss)^(-mu)*(gx/beta_s)); 0.001; 0.001; 0.001]; % initial coefficients for labor parametrization

criter  	= 1e-6;            				        % Convergence criterion
update  	= .65;            				        % Updating rate 

%% Quadrature nodes, weight
n_nodes = quad.n_nodes;
epsi_nodes = quad.epsi_nodes;
weight_nodes = quad.weight_nodes;


%% Loop
iteration = 0; % initializing iteration count
dif = Inf;  % norm(zeta - coef)

while (dif > criter)
    up_bound  = kss * 1.5;     % Upper bound
    low_bound = kss * 0.5;         % Lower bound think this process is widening the bounds as iteration goes on. Thus, limiting the interval of kss at the begining.
for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) * log(kp(t)) +  coef(3) * log(theta(t)) + coef(4) * log(z(t))))^(-1/mu);
    n(t)=n(t)*(n(t)>=N(1))*(n(t)<=N(100))+N(1)*(n(t)<N(1))+N(100)*(n(t)>N(100)); %making n(t) not go over the bounds 
    % calculate u(t) using FOC
    u(t)=((theta(t)*kp(t)^(alpha-1)*n(t)^(1-alpha) - z(t)*phi1)/phi2) +1;
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) / ((1-alpha)*theta(t)*u(t)*kp(t)^(alpha)*n(t)^(-alpha)))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (dss + z(t)* phi1 * (u(t)-1) + phi2/2 * (u(t)-1)^2))*kp(t) + theta(t)*u(t)*kp(t)^(alpha)*n(t)^(1-alpha) - c(t));

    %calculate u(t) counterpart
    if kp(t+1) > up_bound  
          kp(t+1) = up_bound; 
   elseif kp(t+1) < low_bound
          kp(t+1) = low_bound;
    end

end



%% Expectation part calculation
% z(t+1) 
z_p = z.^(rhoz) * exp(sigmaz * epsi_nodes(:,2))';
% theta(t+1)
theta_p = theta.^(rhot) * exp(sigmat * epsi_nodes(:,1))';
%k(t) and k(t+1)
k = repmat(kp(1:end-1),1,n_nodes);
kprime = repmat(kp(2:end),1,n_nodes);
%n(t+1) - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* log(kprime) +  coef(3) .* log(theta_p) + coef(4).* log(z_p))).^(-1/mu));
%u(t+1) - FOC
u_p = (((theta_p.*kprime.^(alpha-1).*n_p.^(1-alpha) - z_p.*phi1)./phi2) + 1);
%c(t+1) - FOC
c_p = ((B.*(1-n_p).^(-mu) ./ ((1-alpha)*theta_p.*u_p.*kprime.^(alpha).*n_p.^(-alpha))).^(-1/gamma));
% Expectation part
e = (1-n_p).^(-mu) ./ (theta.* u.* k.^(alpha).* n.^(-alpha)) .* (theta_p.*u_p.*kprime.^(alpha).*n_p.^(-alpha)) .*...
    ((1- (dss + z_p.* phi1 .* (u_p - 1) + phi2/2 .* (u_p - 1).^2)) + alpha .* theta_p .* u_p .* kprime.^(alpha-1) .* n_p.^(1-alpha)) ;
e = e(1:end-1,:) *weight_nodes;

X = [ones(T-1,1), log(kp(1:end-2)), log(theta(1:end-1)), log(z(1:end-1))];
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


%% Burnout the first 500 periods
N_burnout = 501;
c = c(N_burnout:end,1);
k = k(N_burnout:end,1);
u = u(N_burnout:end,1);
n = n(N_burnout:end,1);
z = z(N_burnout:end,1);
theta = theta(N_burnout:end,1);


c_p = c_p(N_burnout:end,1);
kprime = kprime(N_burnout:end,1);
u_p = u_p(N_burnout:end,1);
n_p = n_p(N_burnout:end,1);
z_p = z_p(N_burnout:end,1);
theta_p = theta_p(N_burnout:end,1);

%% compute residuals - Budget constraint

RBC = ((1- (param.dss + z.* param.phi1 .* (u -1) .* (param.phi2/2) .* (u-1).^2)) .* k + theta .* u .* k.^(param.alpha) .* n.^(1-param.alpha)) ./ (c + param.gx .* kprime) -  1;
REE =  param.beta_s .* c_p.^(-param.gamma) .* (1- (param.dss + z_p.* param.phi1 .* (u_p -1) + param.phi2/2 .* (u_p-1).^2)) + param.alpha .* theta_p .* u_p...
    .* kprime.^(param.alpha-1) .* n_p.^(1-param.alpha) ./ (param.gx .* c.^(-param.gamma)) - 1;

max_RBC = max(log10(abs(RBC)));
mean_RBC = mean(log10(abs(RBC)));

max_REE = max(log10(abs(REE)));
mean_REE = mean(log10(abs(REE)));
