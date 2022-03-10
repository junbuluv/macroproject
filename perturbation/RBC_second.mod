var c n u k z theta y i;
varexo epsiz epsit;

parameters beta alpha dss B gx gamma rhot rhoz sigmat sigmaz mu phi1 phi2 ;
@#include "params.txt"


model;
    [name='Intertemporal Euler Equation']
    gx*c^(-gamma) = beta * gx^(1-gamma) *c(+1)^(-gamma) * ((1- (dss + exp(z(+1))* phi1 * (u(+1)-1) + phi2/2 * (u(+1)-1)^2)) + alpha*exp(theta(+1))*u(+1)*k^(alpha-1)*n(+1)^(1-alpha));
    

    [name='Intratemporal Euler Equation']
    B*(1-n)^(-mu)  = c^(-gamma)*((1-alpha)*exp(theta)*u^(alpha) *k(-1)^(alpha)*n^(-alpha));
    

    [name='Capital utilization']
    (exp(z)*phi1 + phi2 * (u-1))  =  alpha* exp(theta)* u^(alpha-1) * k(-1)^(alpha-1) * n^(1-alpha);
    

    [name='Budget constraint']
    c + gx*k = (1- (dss + exp(z)*phi1 * (u - 1) + phi2/2 * (u-1)^2))*k(-1) + exp(theta)* u^(alpha) * k(-1)^(alpha) * n^(1-alpha);
    

    [name='Production shock']
    theta = rhot*theta(-1)+epsit;
    

    [name='Utilization shock']
    z = rhoz*z(-1)+epsiz;

    [name='Production']
    y = exp(theta) * u * k(-1)^(alpha) * n^(1-alpha);

    [name='Investment']
    i = y - c;



end;

initval;
@#include "steady_state.txt"
end;

steady;
check;

shocks;
var epsit; stderr sigmat;
var epsiz; stderr sigmaz;
end;

stoch_simul(order=2, irf= 60, hp_filter= 1600, periods = 2000, drop = 200);
