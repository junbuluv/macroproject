clear all
clc


CPU = cputime;
%% parameter setting
T = 10000; % simulation length
alpha = 1/3; % capital share to output
nss = 1/3; % labor
k2y = 10; % capital to output ratio
i2y = 0.2133; % investment to output ratio
c2y = 1-i2y; % consumption to output ratio
gamma = 1; % consumption risk aversion
mu = 5; % leisure risk aversion
gx = 1.0029; % Labor augmenting rate

sigmaz   = 0.012; % capital utilization shock standard deviation of disturbance
rhoz = 0.9; % capital utilization shock persistence
sigmat = 0.008; % production shock standard deviation of disturbance
rhot = 0.85; % production shock persistence

dss = 1-gx + (1- c2y) / k2y; % steady state depreciation rate (u = 1)
beta_s = gx / ((1-dss) + alpha * k2y^(-1));
phi1 = gx/beta_s - (1-dss); % capital utilization parameter
phi2 = 2; % capital utilization adjustment parameter
B = (1-alpha) * c2y^(-gamma) * k2y^((alpha*(1-gamma))/(1-alpha)) * nss^(-gamma) * (1-nss)^(mu); % leisure utility parameter
beta = beta_s / gx^(1-gamma);

%% parameter structure
param = struct("beta", beta, "alpha" ,alpha, "mu", mu, "dss",dss, "gx",gx, "phi1",phi1, "phi2",phi2 , "gamma",gamma, "B", B, "rhoz",rhoz, "rhot",rhot,"sigmat",sigmat,"sigmaz",sigmaz);

init = [0.1,0.1,0.33,1];
options = optimset('Display','iter','MaxFunEvals',1000000,'TolFun',1e-8,'MaxIter',10000);
[fsolve_x,fval]=fsolve(@(x) steadycheck(x,param),init,options);

%precomputed steady state values
fsolve_x;


steady_state_names=["c";"k";"n";"u"];
steady_state_vals=[fsolve_x(1); fsolve_x(2); fsolve_x(3); fsolve_x(4)];
writematrix(strcat(strcat(steady_state_names,"="),...
    strcat(string(steady_state_vals),";")),"steady_state.txt");

params_names=["beta";"alpha";"mu";"dss";"gx";"phi1";"phi2";"gamma";"rhot";"rhoz";"sigmat";"sigmaz"; "B"];
params_vals=[param.beta;param.alpha;param.mu;param.dss;param.gx;param.phi1;param.phi2;param.gamma;param.rhot;param.rhoz;param.sigmat;param.sigmaz;param.B];
writematrix(strcat(strcat(params_names,"="),...
    strcat(string(params_vals),";")),"params.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dynare RBC_second.mod


%% get policy function from Dynare


idx_c = 1;
idx_n = 2;
idx_u = 3;
idx_k = 4;
idx_z = 5;
idx_theta = 6;

[A,B,C,D,E,F,G,H,I,J] = get_second_order(idx_c,idx_k,idx_n,idx_u,idx_z,idx_theta,oo_);

% move into a function
% number of nodes we want for each innovation
N_nodes_individual=5;
% number of innovations
N_innov=2;
% variance covariance matrix
vcv_mat=1.0;
[N_nodes_total,qnodes,qweights] = GH_Quadrature(N_nodes_individual,N_innov,vcv_mat);

mat = struct("A",A,"B",B,"C",C,"D",D,"E",E,"F",F,"G",G,"H",H,"I",I,"J",J);

N=10000;
N_states=3;
N_controls=3;
[U,V,Uss,Vss,UE,VE]= simulate_U_V(mat,N,idx_c,idx_k,idx_n,idx_z,idx_theta,idx_u,N_states,N_controls,N_nodes_total,qnodes,oo_,sigmaz,sigmat);


idx_c_v=1;
idx_n_v=2;
idx_u_v=3;
idx_k_u=1;
idx_z_u=2;
idx_theta_u=3;


% Construct residuals
sim_len=length(U(1,:));
E=zeros(sim_len,N_nodes_total);

%To avoid complex number issue apply weight first then
%calculate residual

%%working in progress
%%%%%%%

for j=2:sim_len
    for m=1:N_nodes_total
       E(j,m)=((VE(idx_c_v,j,m).^(-param.gamma).*(((1.0- (param.dss + exp(UE(idx_z_u,j,m)) * param.phi1 * (VE(idx_u_v,j,m) -1) + param.phi2/2 * (VE(idx_u_v,j,m) -1)^2)) +  ...
           param.alpha * exp(UE(idx_theta_u,j,m)).*(U(idx_k_u,j)).^(param.alpha-1) .* VE(idx_u_v,j,m)^(alpha) .* VE(idx_n_v,j,m))^(1-alpha))));
       
    end
end


% apply weights
Eq=E*qweights;

% compute errors
REE=((param.beta * param.gx^(1-param.gamma) .* Eq(2:end,:)) .*(param.gx * (V(idx_c_v,1:end-1).^(-param.gamma))' ).^(-1)) - 1;
RBC = ((V(idx_c_v,1:end-1) + param.gx.*U(idx_k_u,2:end)) ./ ( (1- (param.dss + exp(U(idx_z_u,1:end-1)) .* param.phi1 .* (V(idx_u_v,1:end-1) -1) + param.phi2 .* (V(idx_u_v,1:end-1)-1).^2) ).* U(idx_k_u,1:end-1) ...
      +  exp(U(idx_theta_u,1:end-1)) .* V(idx_u_v,1:end-1) .* U(idx_k_u,1:end-1).^(param.alpha) .* V(idx_n_v,1:end-1).^(1-param.alpha))) -1 ;
RMUL = (param.B .* (1-V(idx_n_v,:)).^(-param.mu)) ./ (V(idx_c_v,:).^(-param.gamma).*  ( (1-param.alpha) .* exp(U(idx_theta_u,:)) .* V(idx_u_v,:).^(param.alpha) .* U(idx_k_u,:).^(param.alpha) .* V(idx_n_v,:).^(-param.alpha) ) ) -1;


max_ee_res=log10(max(abs(REE)));
mean_ee_res=log10(mean(abs(REE)));
max_bc_res=log10(max(abs(RBC)));
mean_bc_res=log10(mean(abs(RBC)));
max_mul_res=log10(max(abs(RMUL)));
mean_mul_res=log10(mean(abs(RMUL)));


fprintf("Max Euler equation residual is %.3f \n",max_ee_res);
fprintf("Mean Euler equation residual is %.3f \n",mean_ee_res);
fprintf("Max Budget Constraint residual is %.3f \n",max_bc_res);
fprintf("Mean Budget Constraint residual is %.3f \n",mean_bc_res);
fprintf("Max Intra-Euler Equation residual is %.3f \n",max_mul_res);
fprintf("Mean Intra-Euler Equation residual is %.3f \n",mean_mul_res);

