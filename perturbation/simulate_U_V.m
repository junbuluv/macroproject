function [U,V,Uss,Vss,UE,VE]=simulate_U_V(mat,N,idx_c,idx_k,idx_n,idx_z,idx_theta,idx_u,N_states,N_controls,N_nodes_total,q_nodes,oo_,sigmaz,sigmat)



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
    rng('default')
    epsi=[randn(1,N+burn_in)*sigmaz;randn(1,N+burn_in)*sigmat];
    
    U=zeros(N_states,N+burn_in);
    V=zeros(N_controls,N+burn_in);
    
    % storage for prime variables need to compute expectations
    UE=zeros(N_states,N+burn_in,N_nodes_total);
    VE=zeros(N_controls,N+burn_in,N_nodes_total);
    
    U(:,1)=B*epsi(:,1) + 0.5*D*(kron(epsi(:,1),epsi(:,1)));
    V(:,1)=G*epsi(:,1) + 0.5*I*(kron(epsi(:,1),epsi(:,1)));
    
    for j=2:N+burn_in
        U(:,j)=A*U(:,j-1)+B*epsi(:,j) + 0.5*C*(kron(U(:,j-1),U(:,j-1))) + 0.5*D*(kron(epsi(:,j),epsi(:,j))) + E*(kron(U(:,j-1),epsi(:,j))); 
        V(:,j)=F*U(:,j-1)+G*epsi(:,j) + 0.5*H*(kron(U(:,j-1),U(:,j-1))) + 0.5*I*(kron(epsi(:,j),epsi(:,j))) + J*(kron(U(:,j-1),epsi(:,j)));
        %%%
        for m=1:N_nodes_total
            UE(:,j,m)=A*U(:,j)+B*q_nodes(m,:)' + 0.5*C*(kron(U(:,j),U(:,j))) + 0.5*D*kron(q_nodes(m,:)',q_nodes(m,:)') + E*(kron(U(:,j),q_nodes(m,:)'));
            VE(:,j,m)=F*U(:,j)+G*q_nodes(m,:)' + 0.5*H*(kron(U(:,j),U(:,j))) + 0.5*I*kron(q_nodes(m,:)',q_nodes(m,:)') + J*(kron(U(:,j),q_nodes(m,:)'));
        end
    end

    Uss=[oo_.dr.ys(idx_k); oo_.dr.ys(idx_z); oo_.dr.ys(idx_theta)];
    Vss=[oo_.dr.ys(idx_c); oo_.dr.ys(idx_n); oo_.dr.ys(idx_u)];

    Uss2 = [oo_.dr.ghs2(idx_k); oo_.dr.ghs2(idx_z); oo_.dr.ghs2(idx_theta)];
    Vss2 = [oo_.dr.ghs2(idx_c); oo_.dr.ghs2(idx_n); oo_.dr.ghs2(idx_u)];

    U=U+Uss*ones(1,N+burn_in) + Uss2*ones(1,N+burn_in);
    V=V+Vss*ones(1,N+burn_in) + Vss2*ones(1,N+burn_in);
    
    UE=UE+Uss*ones(1,N+burn_in) + Uss2*ones(1,N+burn_in);
    VE=VE+Vss*ones(1,N+burn_in) + Vss2*ones(1,N+burn_in);
    
     U=U(:,burn_in+1:end);
     V=V(:,burn_in+1:end);
     UE=UE(:,burn_in+1:end,:);
     VE=VE(:,burn_in+1:end,:);
end