function [err_z,err_p] = errors_exact_HL_k_3_first(p,t,ed,t_e,basis_p1,basis_rt1,z_h,z_vec_r,z_vec_th,z_vec_z,p_h,p_exact,n)
% ERRORS_EXACT_HL_K_3_FIRST - Calculate the errors of our solution x
% compared to the exact solution u.
%   Hodge Laplacian k = 3 case, first order
%
% Syntax:
%     [err_z,err_p] = 
%         errors_exact_HL_k_3_first(p,t,ed,t_e,basis_vertices,basis_rt1,z_h,z_vec_r,z_vec_th,z_vec_z,p_h,p_exact,n)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_e - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     z_h - approximated solution of z
%     z_vec_r - exact solution z vector r component
%     z_vec_r - exact solution z vector theta component
%     z_vec_z - exact solution z vector z component
%     p_h - approximated solution of p
%     p_exact - exact solution p
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%    err_z - L2 error for z approximation
%    err_p - L2 error for p approximation
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
[edges,~] = size(ed);

integral_z = 0;
integral_p = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates
        coordinates(N,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    basis_rt1_r = @(r,z,i) basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
        + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z;
    basis_rt1_th = @(r,z,i) (1./n).*(basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
        + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z);
    basis_rt1_z = @(r,z,i) basis_rt1(4,i,T).*r + basis_rt1(5,i,T).*z ...
        + basis_rt1(6,i,T) + basis_rt1(7,i,T).*r.*z + basis_rt1(8,i,T).*z.*z;

    % basis_p1_r = 0 , basis_p1_z = 0
    basis_p1_th = @(r,z,i) (1./n).*(basis_p1(1,i,T).*r.*r ... 
        + basis_p1(2,i,T).*r.*z + basis_p1(3,i,T).*r);  
    
    approx_r = @(r,z) z_h(nodes+t_e(T,1)).*basis_rt1_r(r,z,1) + z_h(nodes+edges+t_e(T,1)).*basis_rt1_r(r,z,2)...
        + z_h(nodes+t_e(T,2)).*basis_rt1_r(r,z,3) + z_h(nodes+edges+t_e(T,2)).*basis_rt1_r(r,z,4)...
        + z_h(nodes+t_e(T,3)).*basis_rt1_r(r,z,5) + z_h(nodes+edges+t_e(T,3)).*basis_rt1_r(r,z,6)...
        + z_h(nodes+2*edges+T).*basis_rt1_r(r,z,7)...
        + z_h(nodes+2*edges+triangles+T).*basis_rt1_r(r,z,8);
    
    approx_th = @(r,z) z_h(nodes+t_e(T,1)).*basis_rt1_th(r,z,1) + z_h(nodes+edges+t_e(T,1)).*basis_rt1_th(r,z,2)...
        + z_h(nodes+t_e(T,2)).*basis_rt1_th(r,z,3) + z_h(nodes+edges+t_e(T,2)).*basis_rt1_th(r,z,4)...
        + z_h(nodes+t_e(T,3)).*basis_rt1_th(r,z,5) + z_h(nodes+edges+t_e(T,3)).*basis_rt1_th(r,z,6)...
        + z_h(nodes+2*edges+T).*basis_rt1_th(r,z,7)...
        + z_h(nodes+2*edges+triangles+T).*basis_rt1_th(r,z,8)...
        + z_h(t(1,T)).*basis_p1_th(r,z,1)...
        + z_h(t(2,T)).*basis_p1_th(r,z,2)...
        + z_h(t(3,T)).*basis_p1_th(r,z,3); %2*edges+2*triangles+
    
    approx_z = @(r,z) z_h(nodes+t_e(T,1)).*basis_rt1_z(r,z,1)+z_h(nodes+edges+t_e(T,1)).*basis_rt1_z(r,z,2)...
        + z_h(nodes+t_e(T,2)).*basis_rt1_z(r,z,3) + z_h(nodes+edges+t_e(T,2)).*basis_rt1_z(r,z,4)...
        + z_h(nodes+t_e(T,3)).*basis_rt1_z(r,z,5) + z_h(nodes+edges+t_e(T,3)).*basis_rt1_z(r,z,6)...
        + z_h(nodes+2*edges+T).*basis_rt1_z(r,z,7)...
        + z_h(nodes+2*edges+triangles+T).*basis_rt1_z(r,z,8);
       
    % find L2 Error
    integrandz =@(r,z) ((z_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (z_vec_th(r,z) - approx_th(r,z)).^2 ...
        + (z_vec_z(r,z) - approx_z(r,z)).^2).*r;
    
    integral_z = integral_z + Wx'*feval(integrandz,X,Y)*Wy;
    
    % find L2 Error for p
    basis_p1_fn =@(r,z,i) basis_p1(1,i,T).*r ... 
        + basis_p1(2,i,T).*z + basis_p1(3,i,T);
    
    approx = @(r,z) p_h(t(1,T)).*basis_p1_fn(r,z,1)...
        + p_h(t(2,T)).*basis_p1_fn(r,z,2) + p_h(t(3,T)).*basis_p1_fn(r,z,3);
    
    integrandp =@(r,z) ((p_exact(r,z) - approx(r,z)).^2).*r;
    
    integral_p = integral_p + Wx'*feval(integrandp,X,Y)*Wy;        
end

err_z = sqrt(integral_z);
err_p = sqrt(integral_p);
% end