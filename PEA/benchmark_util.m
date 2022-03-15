function [res_table, table1, table2] = benchmark_util(T, param, steady_state, node_number, sim_length,sim_num)


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
sigmat = param.sigmat;
rhot = param.rhot;



%% Upper bound and lower bound
N(1)=(nss)*0.5;  % the lower bound on n(t)
N(100)=(nss)*2;  % the upper bound on n(t)


%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1) = kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
u = zeros(T,1); % u(t)


coef = [log((1-nss)^(-mu)*(gx/beta_s)); 0.001; 0.001]; % initial coefficients for labor parametrization

criter  	= 1e-6;            				        % Convergence criterion
update  	= .6;            				        % Updating rate 

epsi_number = 1;
weight = 1;
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(node_number,epsi_number,weight') ;


%% Loop
iteration = 0; % initializing iteration count
dif = Inf;  % norm(zeta - coef)

while (dif > criter)
    up_bound  = kss * 1.5;     % Upper bound 150% of the steady state kss.
    low_bound = kss * 0.5;         % Lower bound 50% of the steady state kss.
for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) *kp(t) +  coef(3) * theta(t)))^(-1/mu);
    
    n(t)=n(t)*(n(t)>=N(1))*(n(t)<=N(100))+N(1)*(n(t)<N(1))+N(100)*(n(t)>N(100)); %making n(t) not go over the bounds 
    % calculate u(t) using FOC
    u(t) = (alpha * theta(t) * kp(t)^(alpha-1) * n(t)^(1-alpha) / dss )^(1/(phi-alpha));
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) * (1/((1-alpha)*theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(-alpha))))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (dss/phi * u(t)^(phi)))*kp(t) + theta(t)*(u(t)*kp(t))^(alpha)*n(t)^(1-alpha) - c(t));

    %calculate u(t) counterpart
    if kp(t+1) > up_bound  
          kp(t+1) = up_bound; 
   elseif kp(t+1) < low_bound
          kp(t+1) = low_bound;
    end

end



%% Expectation part calculation 
% theta(t+1)
theta_p = repmat(theta.^(rhot),1,n_nodes) .* exp(sigmat * epsi_nodes)';
%k(t) and k(t+1)
k = repmat(kp(1:end-1),1,n_nodes);
kprime = repmat(kp(2:end),1,n_nodes);
%n(t+1) - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* kprime +  coef(3) .* theta_p)).^(-1/mu));
%u(t+1) - FOC
u_p = (alpha .* theta_p .* kprime.^(alpha-1) .* n_p.^(1-alpha) ./ ( dss)).^(1/(phi-alpha));
% Expectation part
e = (1-n_p).^(-mu) .* (theta.* (u.* k).^(alpha).* n.^(-alpha)) ./ (theta_p.*(u_p.*kprime).^(alpha).*n_p.^(-alpha)) .*...
    ((1- (dss/phi .* u_p.^(phi))) + alpha .* theta_p .* u_p.^(alpha) .* kprime.^(alpha-1) .* n_p.^(1-alpha)) ;
e = e(1:end-1,:) *weight_nodes;

X = [ones(T-1,1), kp(1:end-2,:), theta(1:end-1,:)];
zeta = nlinfit(X,e,'object',coef);
dif = norm(coef - zeta);

if rem(iteration,100) == 0
    dif
end
if dif > criter
    coef = update*zeta + (1-update)*coef;
%else
   % fprintf("Optimal Coefficient for parametrized function is %.3f \n", coef);
end
iteration = iteration+1;

end



Z = sim_length;
theta = Shocks(Z,sigmat,rhot);
k_sim = zeros(Z+1,1);
k_sim(1,1) = kss;
n_sim = zeros(Z,1);
c_sim = zeros(Z,1);
u_sim = zeros(Z,1);
for t = 1:Z
   n_sim(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) *k_sim(t) +  coef(3) * theta(t)))^(-1/mu);
   u_sim(t) = (alpha * theta(t) * k_sim(t)^(alpha-1) * n_sim(t)^(1-alpha) / dss )^(1/(phi-alpha));
   c_sim(t) = (B*(1-n_sim(t))^(-mu) * (1/((1-alpha)*theta(t)*(u_sim(t)*k_sim(t))^(alpha)*n_sim(t)^(-alpha))))^(-1/gamma);
   k_sim(t+1) = (1/gx)*((1- (dss/phi * u_sim(t)^(phi)))*k_sim(t) + theta(t)* k_sim(t)^(alpha)*n_sim(t)^(1-alpha) - c_sim(t));

end

theta_p = repmat(theta.^(rhot),1,n_nodes) .* exp(sigmat * epsi_nodes)';
%k(t) and k(t+1)
k = repmat(k_sim(1:end-1),1,n_nodes);
kprime = repmat(k_sim(2:end),1,n_nodes);
%n(t+1) - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* kprime +  coef(3) .* theta_p)).^(-1/mu));
%u(t+1)
u_p = (alpha .* theta_p .* kprime.^(alpha-1) .* n_p.^(1-alpha) ./ (dss)).^(1/(phi-alpha));
%c(t+1) - FOC
c_p = ((B.*(1-n_p).^(-mu) ./ ((1-alpha)*theta_p.*(kprime).^(alpha).*n_p.^(-alpha))).^(-1/gamma));


%% compute residuals 
uniform_weight = ones(n_nodes,1)./n_nodes;
RBC = ((1- param.dss) .* k + theta .* u_sim.^param.alpha .* k.^(param.alpha) .* n_sim.^(1-param.alpha)) ./ (c_sim + param.gx .* kprime) * uniform_weight  -  1;
REE =  (param.beta_s .* c_p.^(-param.gamma) .* ((1- (param.dss/param.phi .* u_p.^param.phi)) + param.alpha .* theta_p .* u_p.^(param.alpha) .* kprime.^(param.alpha-1) .* n_p.^(1-param.alpha)) ./ (param.gx .* c_sim.^(-param.gamma)))*weight_nodes  - 1;
RMUL = ((c_sim.^(-param.gamma) .* (1-param.alpha).*theta.*(u_sim.*k).^(param.alpha).*n_sim.^(-param.alpha))./ (param.B.*(1-n_sim).^(-param.mu))) * uniform_weight  -1;


max_RBC = max(log10(abs(RBC)));
mean_RBC = mean(log10(abs(RBC)));
max_REE = max(log10(abs(REE)));
mean_REE = mean(log10(abs(REE)));
max_RMUL = max(log10(abs(RMUL)));
mean_RMUL = mean(log10(abs(RMUL)));

res_table = table(max_REE,mean_REE,max_RBC,mean_RBC,max_RMUL,mean_RMUL);

k_simu = zeros(Z+1,sim_num);
k_simu(1,:) = kss;
n_simu = zeros(Z,sim_num);
c_simu = zeros(Z,sim_num);
u_simu = zeros(Z,sim_num);
y_simu = zeros(Z,sim_num);
w_simu = zeros(Z,sim_num);
i_simu = zeros(Z,sim_num);


for j= 1:sim_num
theta = Shocks(Z,sigmat,rhot);
for t = 1:Z
   n_simu(t,j) = 1- (beta_s/gx * exp(coef(1) + coef(2) *k_simu(t,j) +  coef(3) * theta(t)))^(-1/mu);
   u_simu(t,j) = (alpha * theta(t) * k_simu(t)^(alpha-1) * n_simu(t)^(1-alpha) / dss )^(1/(phi-alpha));
   c_simu(t,j) = (B*(1-n_simu(t,j))^(-mu) * (1/((1-alpha)*theta(t)*(u_simu(t,j)*k_simu(t,j))^(alpha)*n_simu(t)^(-alpha))))^(-1/gamma);
   k_simu(t+1,j) = (1/gx)*((1- (dss/phi * u_simu(t)^(phi)))*k_simu(t,j) + theta(t)*k_simu(t,j)^(alpha)*n_simu(t,j)^(1-alpha) - c_simu(t,j));
   y_simu(t,j) = theta(t)*(u_simu(t)*k_simu(t))^(alpha)*n_simu(t)^(1-alpha);
   w_simu(t,j) = (1-alpha)*theta(t)*(u_simu(t,j)*k_simu(t,j))^(alpha)*n_simu(t,j)^(-alpha);
   i_simu(t,j) = y_simu(t,j) - c_simu(t,j);
end
end





k_sim = k_simu;
n_sim = n_simu;
c_sim = c_simu;
y_sim = y_simu;
w_sim = w_simu;
i_sim = i_simu;
u_sim = u_simu;



t= (1:1:Z);
t_k= (1:1:Z+1);
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

Interval_model_up = Model+2*[std_pk, std_pi, std(mean(n_sim)), std_sdc, std_sdn, std_sdk, std_sdi, std_sdy, std_sdw, std_sdu, std(diag(corr(y_c(:,t), c_c(:,t)))), std(diag(corr(y_c(:,t), n_c(:,t))))...
    std(diag(corr(y_c(:,t), k_c(:,t)))), std(diag(corr(y_c(:,t), i_c(:,t)))), std(diag(corr(y_c(:,t), w_c(:,t)))), std(diag(corr(y_c(:,t), u_c(:,t))))]';

Interval_model_low = Model-2*[std_pk, std_pi, std(mean(n_sim)), std_sdc, std_sdn, std_sdk, std_sdi, std_sdy, std_sdw, std_sdu, std(diag(corr(y_c(:,t), c_c(:,t)))), std(diag(corr(y_c(:,t), n_c(:,t))))...
    std(diag(corr(y_c(:,t), k_c(:,t)))), std(diag(corr(y_c(:,t), i_c(:,t)))), std(diag(corr(y_c(:,t), w_c(:,t)))), std(diag(corr(y_c(:,t), u_c(:,t))))]';


table1 = table(US_data,Model,Model_std,Interval_model_up, Interval_model_low);
table2 = table(Mean_data,Interval_high,Interval_low);


