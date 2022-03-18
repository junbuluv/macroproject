function [U,V,theta,z]=simulation_per(mat,N,idx_c,idx_k,idx_n,idx_z,idx_theta,idx_u,N_states,N_controls,oo_,sigmaz,sigmat,param,idx)



A = mat.A; 
B = mat.B; 
C = mat.C; 
D = mat.D; 
E = mat.E;
F = mat.F;
G = mat.G;
H = mat.H;
I = mat.I;
J = mat.J;


    burn_in=500;

    epsi=[randn(1,N+burn_in)*sigmaz;randn(1,N+burn_in)*sigmat];
    
    U=zeros(N_states,N+burn_in);
    V=zeros(N_controls,N+burn_in);
    
    % storage for prime variables need to compute expectations
    
    
    U(:,1)=B*epsi(:,1) + 0.5*D*(kron(epsi(:,1),epsi(:,1)));
    V(:,1)=G*epsi(:,1) + 0.5*I*(kron(epsi(:,1),epsi(:,1)));
    
    for j=2:N+burn_in
        U(:,j)=A*U(:,j-1)+B*epsi(:,j) + 0.5*C*(kron(U(:,j-1),U(:,j-1))) + 0.5*D*(kron(epsi(:,j),epsi(:,j))) + E*(kron(U(:,j-1),epsi(:,j))); 
        V(:,j)=F*U(:,j-1)+G*epsi(:,j) + 0.5*H*(kron(U(:,j-1),U(:,j-1))) + 0.5*I*(kron(epsi(:,j),epsi(:,j))) + J*(kron(U(:,j-1),epsi(:,j)));
        %%%
    end

    Uss=[oo_.dr.ys(idx_k); oo_.dr.ys(idx_z); oo_.dr.ys(idx_theta)];
    Vss=[oo_.dr.ys(idx_c); oo_.dr.ys(idx_n); oo_.dr.ys(idx_u)];

    Uss2 = [oo_.dr.ghs2(idx_k); oo_.dr.ghs2(idx_z); oo_.dr.ghs2(idx_theta)];
    Vss2 = [oo_.dr.ghs2(idx_c); oo_.dr.ghs2(idx_n); oo_.dr.ghs2(idx_u)];

    U=U+Uss*ones(1,N+burn_in) + Uss2*ones(1,N+burn_in);
    V=V+Vss*ones(1,N+burn_in) + Vss2*ones(1,N+burn_in);
    
     U=U(:,burn_in+1:end);
     V=V(:,burn_in+1:end);
     theta = epsi(1,burn_in+1:end)';
     z = epsi(2,burn_in+1:end)';
end