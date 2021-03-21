function [basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h] = weighted_HL_k_2_first(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_2_FIRST Hodge Laplacian k = 1 firt order Finite Element Method.
%   This program is set up to give approximations of the unknown solutions 
%   u and s.
%   Hodge Laplacian k = 2 case, first order
%   {zeta_j}j=1->N+Ne is the basis for Bh1
%   {psi_i}i=1->Ne+Nt is the basis for Ch1
%   (Bh1 is the weighted fourier ND1 and P2 space)
%   (Ch1 is the weighted fourier RT1 and P1 space)
%   Solve for (s,u) in (Bh x Ch) s.t.
%       (s , w)_r - (curl_rz^n(w) , u)_r = 0
%       (curl_rz^n(s) , v)_r + (div_rz^n(u) , div_rz^n(v))_r = (f , v)_r
%           for all w in Bh, v in Ch
%
% Syntax:
%     [basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h]
%       = weighted_HL_k_2_first(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n)
%     f_vec_r - given function r component
%     f_vec_th - given function theta component
%     f_vec_z - given function z component
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T)
%         represents the pieceiwise basis function for the ith edge in
%         triangle T.
%     basis_p2 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_p2(:,i,T)
%         represents the pieceiwise basis function for the ith node or 
%         midpoint in triangle T.
%     basis_nd1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_nd1(:,i,T)
%         represents the ith pieceiwise basis function in triangle T.
%     u_h - approximated solution for u
%     s_h - approximated solution for s
%
% Usage Exampled:
%    addpath ../helper_functions data
%    n = 1;
%    [~,~,~,~,~,~,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h] = weighted_HL_k_2_first(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n);
% Dependencies:
%    ../new_ele[mesh].mat
%    basis_functions_HL_k_2_first.m
%    create_B_HL_k_2_first.m
%    create_F_HL_k_2_first.m
%    create_S_HL_k_2_first.m
%    mass_matrix_HL_k_2_first.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../');
addpath('../edge_resources/');

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);

% To ensure we refine every triangle the same
[~,num_node]=size(p);
it=zeros(1,num_node);
for i=1:num_node
    it(i)=i;
end   

for i = 2:mesh
    fprintf('%d\t',i);
    % Refine mesh to next level
    [p,e,t]=refinemesh(g,p,e,t,it,'regular');    
end

[~,triangles]=size(t);
ele=t(1:3,1:triangles);
ele=ele';
node=p';
tr=triangulation(ele,node);
ed = edges(tr);
%load(['../new_ele',num2str(mesh),'.mat']);
%t_ed = new_ele
load(['t_ed_',num2str(mesh),'.mat']);
t_ed = t_ed';

% Find the midpoints for P2 nodal points
[p2,t2] = find_midpoints(p,t);

[basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h]...
    = solve(p,t,p2,t2,ele,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);

% end main
end

% subfunction
function [basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h] = solve(p,t,p2,t2,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n)
    [basis_p1,basis_rt1,basis_p2,basis_nd1] = basis_functions_HL_k_2_first(p,t,p2,t2,ele,ed,new_ele);
    M = mass_matrix_HL_k_2_first(p,t,p2,t2,ed,new_ele,basis_p2,basis_nd1,n);
    B = create_B_HL_k_2_first(p,t,p2,t2,ed,new_ele,basis_p1,basis_rt1,basis_p2,basis_nd1,n);
    S = create_S_HL_k_2_first(p,t,ed,new_ele,basis_p1,basis_rt1);
    F = create_F_HL_k_2_first(p,t,ed,new_ele,basis_p1,basis_rt1,f_vec_r,f_vec_th,f_vec_z,n);
    
    Bt = B.';
    MB = (M\B);
    u_h = (Bt*MB + S)\F;
    s_h = M\(B*u_h);
end