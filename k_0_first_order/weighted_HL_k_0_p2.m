function [basis,u_h] = weighted_HL_k_0_p2(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_0_P2 Weighted Hodge Laplacian with k = 0. 
%   This program is set up to give an approximation of an unknown solution 
%   u.
%
% Syntax:
%     [basis,u_h] = weighted_HL_k_0_p2(u,grad_u_r,grad_u_z,gd,sf,ns,mesh,n)
%
% Inputs:
%     f - given function
%     grad_f_r - gradient(f) w.r.t r
%     grad_f_z - gradient(f) w.r.t z
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     basis - a 6x6xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     u_h - approximated solution for u
%
% Usage Exampled:
%    addpath ../../data ../data/
%    [f,grad_f_r,grad_f_z,~,~,~] = get_data7();
%    mesh = 8;
%    n = 1;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [basis,u_h] = weighted_HL_k_0_p2(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n);
% Dependencies:
%    find_midpoints.m
%    basis_functions_weighted_HL_k_0_p2.m
%    create_b_HL_k_0_p2.m
%    stiffness_matrix_weighted_HL_k_0_p2.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../../')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
pdemesh(p,e,t, 'NodeLabels','off', 'ElementLabels','off');

% To ensure we refine every triangle the same
[~,num_node]=size(p);
it=zeros(1,num_node);
for i=1:num_node
    it(i)=i;
end

for i = 2:mesh
    % Refine mesh to next level
    [p,e,t]=refinemesh(g,p,e,t,it,'regular');
    %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');
end
pdemesh(p,e,t, 'NodeLabels','off', 'ElementLabels','off');

% Find the midpoints for P2 nodal points
[p2,t2] = find_midpoints(p,t);

% solve
[basis,u_h] = solve(p,p2,e,t,t2,f,grad_f_r,grad_f_z,n);

% end main
end

% subfunction
function [basis,u_h] = solve(p,p2,e,t,t2,f,grad_f_r,grad_f_z,n)
    basis = basis_functions_weighted_HL_k_0_p2(p,t,p2,t2);
    S = stiffness_matrix_weighted_HL_k_0_p2(p,t,p2,t2,basis,n);
    b = create_b_HL_k_0_p2(p,t,p2,t2,basis,f,grad_f_r,grad_f_z,n);
    u_h = S\b;

    figure();
  %  pdeplot([p,p2],e,t, 'XYData',u_h, 'ZData', u_h, 'Mesh', 'on');
end
