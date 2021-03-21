function [basis_p1,basis_rt1,z_h,p_h] = weighted_HL_k_3_first(f,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_3_FIRST Hodge Laplacian k = 3 First Order Finite Element Method.
%   This program is set up to give approximations of the unknown solutions 
%   z and p.
%   Hodge Laplacian k = 3 case, first order
%   {psi_i}i=1->2Ne+2Nt+N is the basis for Ch1
%   {chi_j}j=1->2Nt is the basis for Dh1
%   (Ch is the weighted fourier Raviart Thomas space)
%   (Dh is the piecewise constant space)
%   Solve for (z,p) in (Ch x Dh) s.t.
%       (z , w)_r - (p , div_rz^n(w))_r = 0
%       (div_rz^n(z) , s)_r = (f , s)_r
%           for all w in Ch, s in Dh
%
% Syntax:
%     [basis_p1,basis_rt1,z_h,p_h] = weighted_HL_k_3_first(f,gd,sf,ns,mesh,n)
%     f - given function
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     basis_vertices - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     z_h - approximated solution for z
%     p_h - approximated solution for p
%
% Usage Exampled:
%    addpath ../helper_functions data
%    n = 1;
%    [~,~,~,~,f] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [basis_p1,basis_rt1,z_h,p_h] = weighted_HL_k_3_first(f,gd,sf,ns,mesh,n);
% Dependencies:
%    ../new_ele[mesh].mat.m
%    basis_functions_HL_k_3_first.m
%    create_B_HL_k_3_first.m
%    create_F_HL_k_3_first.m
%    display_errors.m
%    errors_exact_HL_k_3_first.m
%    mass_matrix_HL_k_3_first.m
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

% solve
[basis_p1,basis_rt1,z_h,p_h] = solve(p,t,ele,ed,t_ed,f,n);

% end main
end

% subfunction
function [basis_p1,basis_rt1,z_h,p_h] = solve(p,t,ele,ed,new_ele,f,n)
    [basis_p1,basis_rt1] = basis_functions_HL_k_3_first(p,t,ele,ed,new_ele);
    M = mass_matrix_HL_k_3_first(p,t,ed,new_ele,basis_p1,basis_rt1,n);
    B = create_B_HL_k_3_first(p,t,ed,new_ele,basis_p1,basis_rt1);
    F = create_F_HL_k_3_first(p,t,basis_p1,f);
    
    Bt = B.';
    MB = (M\B);
    p_h = (Bt*MB)\F;
    z_h = MB*p_h;
end