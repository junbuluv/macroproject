function tet = Shocks(T,sigma,rho)
tet     = zeros(T,1)+1;                 
epsi    = randn(T,1)*sigma;
for t=2:T
   tet(t)=tet(t-1)^rho*exp(epsi(t));
end
