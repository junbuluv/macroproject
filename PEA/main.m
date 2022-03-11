clear all
clc


CPU = cputime;
%% parameter setting
T = 10000; % simulation length
alpha = 1/3; % capital share to output
nss = 1/3; % labor steady state
k2y = 10; % capital to output ratio
i2y = 0.2133; % investment to output ratio
c2y = 1-i2y; % consumption to output ratio
gamma = 1; % consumption risk aversion
mu = 5; % leisure risk aversion
gx = 1.0029; % Labor augmenting rate
sigmat = 0.016; % production shock standard deviation of disturbance
rhot = 0.95; % production shock persistence
uss = 1; % utilization rate in steady state
phi = 1.81; % capital utilization parameter
dss = alpha * k2y^(-1) / uss^(phi); % steady state depreciation rate (u = 1)
beta_s = gx / ((1- dss/phi * uss^(phi)) + alpha * k2y^(-1));
B = c2y^(-gamma) * (1-alpha) * k2y^(alpha*(1-gamma)/(1-alpha)) * uss^(alpha*(1-gamma)/(1-alpha)) * nss^(-gamma) * (1-nss)^(mu) ; % leisure utility parameter

%% steady state value calculation
param = struct("beta_s", beta_s, "alpha" ,alpha, "mu", mu, "dss",dss, "gx",gx, "phi",phi , "gamma",gamma,...
    "B", B, "rhot",rhot, "sigmat",sigmat);
init = [0.1,0.1,1/3,0.8];
options = optimset('Display','iter','MaxFunEvals',1000000,'TolFun',1e-8,'MaxIter',10000);
[ss_val,f_val]=fsolve(@(x) steady_check(x,param),init,options);


%% steady state
css = ss_val(1);
kss = ss_val(2);
nss = ss_val(3);
uss = ss_val(4);
% save steady state
steady_state = struct("css",css,"kss",kss,"nss",nss,"uss",uss);





%% Quadrature nodes
node_number = 5;
epsi_number = 2;
weight = diag(ones(epsi_number,1));
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(node_number,epsi_number,weight') ;
quad = struct("n_nodes",n_nodes,"epsi_nodes",epsi_nodes,"weight_nodes",weight_nodes);
sigmaz = 0.0013;
rhoz = 0.95;


%% PEA fitting get coefficients for n(t)
coef = fitting(sigmaz,rhoz,T,param,steady_state,quad);

%% compute residual
N_burnout = 500;
[mean_REE, max_REE, mean_RBC, max_RBC, mean_RMUL, max_RMUL] = residual(sigmaz,rhoz,coef,T,N_burnout,param,steady_state,quad)



%% fminbnd (fixed persistence)
options_fminbnd = optimset('Display','iter','MaxFunEvals',1000,'TolFun',1e-10,'MaxIter',5000);
fun2 = @(x) residual(x, 0.95 ,coef,T,N_burnout,param ,steady_state,quad);
[val, fval, exitflag, output] = fminbnd(fun2, 0.0000001,0.5,options_fminbnd)


%% simulation
sigmaz = 0.0012;
rhoz = 0.95;
N = 1172; % simulation length
Z = 5000; % repetition length
Burn = 100;

% initialization

c_sim = zeros(N - Burn,Z);
k_sim = zeros(N+1 - Burn,Z);
n_sim = zeros(N - Burn,Z);
u_sim = zeros(N - Burn,Z);
y_sim = zeros(N - Burn,Z);
w_sim = zeros(N - Burn,Z);
i_sim = zeros(N - Burn,Z);
for t = 1:Z
[c_sim(:,t),k_sim(:,t),n_sim(:,t),u_sim(:,t), y_sim(:,t), w_sim(:,t), i_sim(:,t)] = simulation(sigmaz, rhoz, coef , N, Burn, param, steady_state);
end




corr_yc = mean((diag(corr(y_sim,c_sim))));
corr_yn = mean((diag(corr(y_sim,n_sim))));
corr_yi = mean((diag(corr(y_sim,i_sim))));
corr_yw = mean((diag(corr(y_sim,w_sim))));
corr_yk = mean((diag(corr(y_sim,k_sim(1:end-1,:)))));
corr_apl = mean((diag(corr(y_sim./n_sim,n_sim))));


t= (1:1:N-Burn);
t_k= (1:1:N-Burn+1);
gx_t = param.gx.^(t)';
gx_t_k = param.gx.^(t_k)';


k_sim = k_sim .* gx_t_k;
c_sim = c_sim .* gx_t;
y_sim = y_sim .* gx_t;
i_sim = y_sim - c_sim;
w_sim = w_sim .* gx_t;

%% calculate ratio  
% capital to output
pk = k_sim(1:end-1,:) ./ y_sim;
% consumption to output
pc = c_sim ./ y_sim;
% investment to output
pi = i_sim ./ y_sim;

% calculate average of the ratio
pk_ss = mean(mean(pk));
pc_ss = mean(mean(pc));
pi_ss = mean(mean(pi));

%% compute the logarithm of the series
ln_c = log(c_sim);
ln_k = log(k_sim);
ln_n = log(n_sim);
ln_u = log(u_sim);
ln_y = log(y_sim);
ln_i = log(i_sim);
ln_w = log(w_sim);

% Hodrick Prescott filter
t= (1:1:Z);

[c_t(:,t), c_c(:,t)] = hpfilter(ln_c(:,t),1600);
[k_t(:,t), k_c(:,t)] = hpfilter(ln_k(:,t),1600);
[n_t(:,t), n_c(:,t)] = hpfilter(ln_n(:,t),1600);
[u_t(:,t), u_c(:,t)] = hpfilter(ln_u(:,t),1600);
[y_t(:,t), y_c(:,t)] = hpfilter(ln_y(:,t),1600);
[i_t(:,t), i_c(:,t)] = hpfilter(ln_i(:,t),1600);
[w_t(:,t), w_c(:,t)] = hpfilter(ln_w(:,t),1600);




%compute the standard deviations



sd_c(:,t) = sqrt(var(c_c(:,t)));
sd_k(:,t) = sqrt(var(k_c(:,t)));
sd_n(:,t) = sqrt(var(n_c(:,t)));
sd_u(:,t) = sqrt(var(u_c(:,t)));
sd_i(:,t) = sqrt(var(i_c(:,t)));
sd_y(:,t) = sqrt(var(y_c(:,t)));
sd_w(:,t) = sqrt(var(w_c(:,t)));


% compute correlation of each variable with output

corr_yc(:,t) = corr(y_c(:,t), c_c(:,t));
corr_yn(:,t) = corr(y_c(:,t), n_c(:,t));




