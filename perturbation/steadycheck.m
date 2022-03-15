function err=steadycheck(init, mp)

beta= mp.beta;
alpha = mp.alpha;
mu = mp.mu;
dss = mp.dss;
gx = mp.gx;
phi = mp.phi;
B = mp.B;
gamma = mp.gamma;


c = init(1);
k = init(2);
n = init(3);
u = init(4);


err(1) = beta *gx^(-gamma) * ((1- (dss/phi * u^phi)) + alpha*u^(alpha)*k^(alpha-1)*n^(1-alpha)) - 1; %intertemporal foc
err(2) = B*(1-n)^(-mu) -(c^(-gamma) * ((1-alpha) * u^(alpha) *k^(alpha) * n^(-alpha))); %intratemporal foc
err(3) = dss*u^(phi-1) - (alpha * u^(alpha-1)* k^(alpha-1) * n^(1-alpha)); % cap utilization foc
err(4) = c + gx*k - (1- dss/phi * u^(phi))*k - u^(alpha)*k^(alpha)*n^(1-alpha); % budget constraint



end

