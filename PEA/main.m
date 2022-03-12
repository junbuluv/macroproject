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
sim_length = 72;
%% PEA fitting get coefficients for n(t)
coef = fitting(sigmaz,rhoz,T,param,steady_state,quad);

%%  PEA benchmark_uilization with no shock
coef_benchmark = benchmark_util(T,param,steady_state,node_number,sim_length);

%% PEA benchmark_no utilization 
coef_no_u = benchmark_no_u(T,param,steady_state,node_number,sim_length);


%% compute residual
N_burnout = 500;
[mean_REE, max_REE, mean_RBC, max_RBC, mean_RMUL, max_RMUL] = residual(sigmaz,rhoz,coef,T,N_burnout,param,steady_state,quad)



%% fminbnd (fixed persistence)
options_fminbnd = optimset('Display','iter','MaxFunEvals',1000,'TolFun',1e-10,'MaxIter',5000);
fun2 = @(x) sigma_search(x, 0.95,T,N_burnout,param ,steady_state,quad);
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
corr_yw = mean((diag(corr(y_sim,w_sim))));
corr_yk = mean((diag(corr(y_sim,k_sim(1:end-1,:)))));
corr_yu = mean((diag(corr(y_sim,u_sim))));
corr_yi = mean((diag(corr(y_sim,i_sim))));
t= (1:1:N-Burn);
t_k= (1:1:N-Burn+1);
gx_t = param.gx.^(t)';
gx_t_k = param.gx.^(t_k)';


k_sim_gx = k_sim .* gx_t_k;
c_sim_gx = c_sim .* gx_t;
y_sim_gx = y_sim .* gx_t;
i_sim_gx = y_sim_gx - c_sim_gx;
w_sim_gx = w_sim .* gx_t;

%% calculate ratio  
% capital to output
pk = k_sim_gx(1:end-1,:) ./ y_sim_gx;
% consumption to output
pc = c_sim_gx ./ y_sim_gx;
% investment to output
pi = i_sim_gx ./ y_sim_gx;

% calculate average of the ratio
pk_ss = mean(mean(pk));
pc_ss = mean(mean(pc));
pi_ss = mean(mean(pi));
std_pk = std(mean(pk));
std_pc = std(mean(pc));
std_pi = std(mean(pi));
%% compute the logarithm of the series
ln_c = log(c_sim_gx);
ln_k = log(k_sim_gx);
ln_n = log(n_sim);
ln_u = log(u_sim);
ln_y = log(y_sim_gx);
ln_i = log(i_sim_gx);
ln_w = log(w_sim_gx);

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



sd_c = mean(sqrt(var(c_c(:,t))));
sd_k = mean(sqrt(var(k_c(:,t))));
sd_n = mean(sqrt(var(n_c(:,t))));
sd_u = mean(sqrt(var(u_c(:,t))));
sd_i = mean(sqrt(var(i_c(:,t))));
sd_y = mean(sqrt(var(y_c(:,t))));
sd_w = mean(sqrt(var(w_c(:,t))));

std_sdc = std(sqrt(var(c_c(:,t))));
std_sdk = std(sqrt(var(k_c(:,t))));
std_sdn = std(sqrt(var(n_c(:,t))));
std_sdu = std(sqrt(var(u_c(:,t))));
std_sdi = std(sqrt(var(i_c(:,t))));
std_sdy = std(sqrt(var(y_c(:,t))));
std_sdw = std(sqrt(var(w_c(:,t))));
% compute correlation of each variable with output
%truncate one period of kprime
k_c = k_c(2:end,:);
corr_yc_gx = mean(diag(corr(y_c(:,t), c_c(:,t))));
corr_yn_gx = mean(diag(corr(y_c(:,t), n_c(:,t))));
corr_yk_gx = mean(diag(corr(y_c(:,t), k_c(:,t))));
corr_yi_gx = mean(diag(corr(y_c(:,t), i_c(:,t))));
corr_yw_gx = mean(diag(corr(y_c(:,t), w_c(:,t))));
corr_yu_gx = mean(diag(corr(y_c(:,t), u_c(:,t))));

int_corr_yc = [corr_yc_gx - 2*std(diag(corr(y_c(:,t), c_c(:,t)))),corr_yc_gx + 2*std(diag(corr(y_c(:,t), c_c(:,t))))];
int_corr_yn = [corr_yn_gx - 2*std(diag(corr(y_c(:,t), n_c(:,t)))),corr_yn_gx + 2*std(diag(corr(y_c(:,t), n_c(:,t))))];
int_corr_yk = [corr_yk_gx - 2*std(diag(corr(y_c(:,t), k_c(:,t)))),corr_yk_gx + 2*std(diag(corr(y_c(:,t), k_c(:,t))))];
int_corr_yi = [corr_yi_gx - 2*std(diag(corr(y_c(:,t), i_c(:,t)))),corr_yi_gx + 2*std(diag(corr(y_c(:,t), i_c(:,t))))];
int_corr_yw = [corr_yw_gx - 2*std(diag(corr(y_c(:,t), w_c(:,t)))),corr_yw_gx + 2*std(diag(corr(y_c(:,t), w_c(:,t))))];
int_corr_yu = [corr_yu_gx - 2*std(diag(corr(y_c(:,t), u_c(:,t)))),corr_yu_gx + 2*std(diag(corr(y_c(:,t), u_c(:,t))))];


% compute interval for series


int_c =[mean(mean(c_sim))-2*std(mean(c_sim)),mean(mean(c_sim))+2*std(mean(c_sim))];
int_n =[mean(mean(n_sim))-2*std(mean(n_sim)),mean(mean(n_sim))+2*std(mean(n_sim))];
int_k =[mean(mean(k_sim))-2*std(mean(k_sim)),mean(mean(k_sim))+2*std(mean(k_sim))];
int_i =[mean(mean(i_sim))-2*std(mean(i_sim)),mean(mean(i_sim))+2*std(mean(i_sim))];
int_y =[mean(mean(y_sim))-2*std(mean(y_sim)),mean(mean(y_sim))+2*std(mean(y_sim))];
int_w =[mean(mean(w_sim))-2*std(mean(w_sim)),mean(mean(w_sim))+2*std(mean(w_sim))];
int_u =[mean(mean(u_sim))-2*std(mean(u_sim)),mean(mean(u_sim))+2*std(mean(u_sim))];


Model = [pk_ss,pi_ss,mean(mean(n_sim)),sd_c,sd_n,sd_k,sd_i,sd_y,sd_w,sd_u,corr_yc_gx,corr_yn_gx,corr_yk_gx,corr_yi_gx,corr_yu_gx,corr_yw_gx]';
Model_std = [std_pk, std_pi, std(mean(n_sim)), std_sdc, std_sdn, std_sdk, std_sdi, std_sdy, std_sdw, std_sdu, std(diag(corr(y_c(:,t), c_c(:,t)))), std(diag(corr(y_c(:,t), n_c(:,t))))...
    std(diag(corr(y_c(:,t), k_c(:,t)))), std(diag(corr(y_c(:,t), i_c(:,t)))), std(diag(corr(y_c(:,t), w_c(:,t)))), std(diag(corr(y_c(:,t), u_c(:,t))))]';
US_data = [10,0.2133,1/3,0.0085,0.0138,0.0070,0.0426,0.0117,0.0102,0,0.8581,0.8333,-0.21,0.9124,0,0.3483]';
Mean_data = [mean(mean(c_sim)),mean(mean(n_sim)), mean(mean(u_sim)), mean(mean(k_sim)), mean(mean(y_sim)), mean(mean(i_sim)), mean(mean(w_sim))]';
Interval = [int_c; int_n; int_u; int_k; int_y; int_i; int_w]';
Interval_high = Interval(2,:)';
Interval_low = Interval(1,:)';
Table1 = table(Model,US_data);
Table2 = table(Mean_data,Interval_high,Interval_low);

