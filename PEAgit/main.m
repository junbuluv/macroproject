clear all
clc


CPU = cputime;
%% parameter setting
T = 10000; % simulation length
alpha = 0.33; % capital share to output
nss = 0.33; % labor
k2y = 10; % capital to output ratio
i2y = 0.2133; % investment to output ratio
c2y = 1-i2y; % consumption to output ratio
gamma = 1; % consumption risk aversion
mu = 5; % leisure risk aversion
gx = 1.0029; % Labor augmenting rate

sigmaz   = 0.0072; % capital utilization shock standard deviation of disturbance
rhoz = 0.95; % capital utilization shock persistence
sigmat = 0.0096; % production shock standard deviation of disturbance
rhot = 0.95; % production shock persistence

dss = 1-gx + (1- c2y) / k2y; % steady state depreciation rate (u = 1)
beta_s = gx / ((1-dss) + alpha * k2y^(-1));
phi1 = gx/beta_s - (1-dss); % capital utilization parameter
phi2 = 100; % capital utilization adjustment parameter
B = (1-alpha) * c2y^(-gamma) * k2y^((alpha*(1-gamma))/(1-alpha)) * nss^(-gamma) * (1-nss)^(mu); % leisure utility parameter

%% steady state value calculation
param = struct("beta_s", beta_s, "alpha" ,alpha, "mu", mu, "dss",dss, "gx",gx, "phi1",phi1, "phi2",phi2 , "gamma",gamma,...
    "B", B, "rhot",rhot, "rhoz",rhoz, "sigmat",sigmat,"sigmaz",sigmaz);
init = [0.1,0.1,0.33,1];
options = optimset('Display','iter','MaxFunEvals',1000000,'TolFun',1e-8,'MaxIter',10000);
[ss_val,fval]=fsolve(@(x) steadyPEA(x,param),init,options);


%% steady state
ss_val

css = ss_val(1);
kss = ss_val(2);
nss = ss_val(3);
uss = ss_val(4);



%% Upper bound and lower bound
up_n = nss * 0.5;
low_n = nss * 2;

%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1) + kss; % k(t+1)
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
e = zeros(T-1,1); % Parameterized expectation 
zp = zeros(T-1,1); % zprime(t)
cp = zeros(T-1,1); % cprime(t)
u = zeros(T,1); % u(t)


coef = [log((1-nss)^(-mu)*(gx/beta_s)); 0.001; 0.001; 0.001]; % initial coefficients
crate    = 0.007;                                   % Speed of moving bounds
criter  	= 1e-6;            				        % Convergence criterion
update  	= .5;            				        % Updating rate 
maxiter  = 100;                                     % Max number of iterations
options = optimset('MaxFunEvals',1000,'TolFun',1e-7,'MaxIter',10000,'Display','off');

%% Quadrature nodes, weight
node_number = 3;
epsi_number = 2;
weight = diag(ones(epsi_number,1));
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(node_number,epsi_number,weight') ;


%% Loop
iteration = 0; % initializing iteration count
dif = Inf;  % norm(zeta - coef)

%while (dif > criter)||(hit==1)
    up_bound = 0.5 * kss; % set upper bound for k
    low_bound = 1.5 * kss;  % set lower bound for k
    hit = 0;
for t = 1:T
    % Compute the time series
    % calculate n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) * log(kp(t)) +  coef(3) * log(theta(t)) + coef(4) * log(z(t))))^(-1/mu);
    % impose bound for n(t)
    if n(t) > up_n
        n(t) = (alpha * low_n + (1-alpha) * up_n);
    elseif n(t) < low_n
        n(t) = (alpha * up_n + (1-alpha) * low_n);
    end
    %calculate u(t) using FOC
    u(t) = fsolve(@(u) findu(u,param,kp(t),n(t),z(t),theta(t)), uss, options);
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) / ((1-alpha)*theta(t)*u(t)^(alpha)*kp(t)^(alpha)*n(t)^(-alpha)))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = ((1- (dss + phi1 * (u(t)-1) + phi2/2 * (u(t)-1)^2))*kp(t) + theta(t)*u(t)^(alpha)*kp(t)^(alpha)*n(t)^(1-alpha) - c(t)) / gx;
    if kp(t+1) > up_bound
        kp(t+1) = (1-alpha) * up_bound + alpha * low_bound;
    elseif kp(t+1) < low_bound
        kp(t+1) = (1-alpha) * low_bound + alpha * up_bound;
    end

   

end
for t = 1:T
    for m = 1:n_nodes
        zprime(t,m) = z(t)^(rhoz) * exp(sigmaz * epsi_nodes(m,2));
        thetaprime(t,m) = theta(t)^(rhot) * exp(sigmat * epsi_nodes(m,1));
        nprime(t,m) = 1- (beta_s/gx * exp(coef(1) + coef(2) * log(kp(t+1)) +  coef(3) * log(thetaprime(t,m)) + coef(4) * log(zprime(t,m))))^(-1/mu);
        uprime(t,m) = fsolve(@(u) findu(u,param,kp(t+1),nprime(t,m),zprime(t,m),thetaprime(t,m)), uss, options);
        cprime(t,m) = (B*(1-nprime(t,m))^(-mu) / ((1-alpha)*thetaprime(t,m)*uprime(t,m)^(alpha)*kp(t+1)^(alpha)*nprime(t,m)^(-alpha)))^(-1/gamma);
        e(t,m) = cprime(t,m)^(-gamma)*((1 - (dss + phi1 * (uprime(t,m)-1) + phi2 * (uprime(t,m)-1)^2)) + alpha * thetaprime(t,m) * uprime(t,m)^(alpha)*kp(t+1)^(alpha-1)*nprime(t,m)^(1-alpha));
    end

end

%% vectorization
%idx = linspace(1,n_nodes*size(z,1),10);
%idx = round(idx);


%zprime = z * epsi_nodes(:,2)';
%zp = vertcat(zprime(:,1),zprime(:,2),zprime(:,3),zprime(:,4),zprime(:,5),zprime(:,6),zprime(:,7),zprime(:,8),zprime(:,9));
%thetaprime = theta * epsi_nodes(:,1)';
%thetap = vertcat(thetaprime(:,1),thetaprime(:,2),thetaprime(:,3),thetaprime(:,4),thetaprime(:,5),thetaprime(:,6),thetaprime(:,7),thetaprime(:,8),thetaprime(:,9));

eprime = e*weight_nodes

X = [ones(T,1), kp(1:end-1), theta(1:end), z(1:end)];
zeta = nlinfit(X,e,'object',coef);


