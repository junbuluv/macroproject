function [A,B,C,D,E,F,G,H,I,J]=get_second_order(idx_c,idx_k,idx_n,idx_u,idx_z,idx_theta,oo_)

    A=[oo_.dr.ghx(oo_.dr.inv_order_var(idx_k),:);
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_theta),:)];

    B=[oo_.dr.ghu(oo_.dr.inv_order_var(idx_k),:);
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_theta),:)];

    C=[oo_.dr.ghxx(oo_.dr.inv_order_var(idx_k),:);
       oo_.dr.ghxx(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghxx(oo_.dr.inv_order_var(idx_theta),:)];
    
    D=[oo_.dr.ghuu(oo_.dr.inv_order_var(idx_k),:);
       oo_.dr.ghuu(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghuu(oo_.dr.inv_order_var(idx_theta),:)];
    
    E=[oo_.dr.ghxu(oo_.dr.inv_order_var(idx_k),:);
       oo_.dr.ghxu(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghxu(oo_.dr.inv_order_var(idx_theta),:)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F=[oo_.dr.ghx(oo_.dr.inv_order_var(idx_c),:);...
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_u),:)];

    G=[oo_.dr.ghu(oo_.dr.inv_order_var(idx_c),:);...
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_u),:)];

    H=[oo_.dr.ghxx(oo_.dr.inv_order_var(idx_c),:);
       oo_.dr.ghxx(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghxx(oo_.dr.inv_order_var(idx_u),:)];

    I=[oo_.dr.ghuu(oo_.dr.inv_order_var(idx_c),:);
       oo_.dr.ghuu(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghuu(oo_.dr.inv_order_var(idx_u),:)];

    J=[oo_.dr.ghxu(oo_.dr.inv_order_var(idx_c),:);
       oo_.dr.ghxu(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghxu(oo_.dr.inv_order_var(idx_u),:)];
end
