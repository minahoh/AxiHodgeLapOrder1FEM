function [err] = errors_exact_ch_1(p,t,ed,t_e,basis_vertices,basis_rt1,x,u_vec_r,u_vec_th,u_vec_z,n)
% ERRORS_EXACT_CH_1 - Calculate the errors of our solution x
% compared to the exact solution u.
%
% Syntax:
%     [err] = 
%         errors_exact_ch_1(p,t,ed,t_ed,basis_edges,basis_triangles,x,u_vec_r,u_vec_th,u_vec_z,n)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_e - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
%     basis_vertices - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     u_vec_r - exact solution vector r component
%     u_vec_th - exact solution vector theta component
%     u_vec_z - exact solution vector z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%    err - L2 error
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
[edges,~] = size(ed);
integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    

    basis_rt1_r = @(r,z,i) basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
        + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z;
    basis_rt1_th = @(r,z,i) (1./n).*(basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
        + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z);
    basis_rt1_z = @(r,z,i) basis_rt1(4,i,T).*r + basis_rt1(5,i,T).*z ...
        + basis_rt1(6,i,T) + basis_rt1(7,i,T).*r.*z + basis_rt1(8,i,T).*z.*z;

    % basis_p1_r = 0 , basis_p1_z = 0
    basis_p1_th = @(r,z,i) (1./n).*(basis_vertices(1,i,T).*r.*r ... 
        + basis_vertices(2,i,T).*r.*z + basis_vertices(3,i,T).*r);  
    
    approx_r = @(r,z) x(nodes+t_e(T,1)).*basis_rt1_r(r,z,1) + x(nodes+edges+t_e(T,1)).*basis_rt1_r(r,z,2)...
        + x(nodes+t_e(T,2)).*basis_rt1_r(r,z,3) + x(nodes+edges+t_e(T,2)).*basis_rt1_r(r,z,4)...
        + x(nodes+t_e(T,3)).*basis_rt1_r(r,z,5) + x(nodes+edges+t_e(T,3)).*basis_rt1_r(r,z,6)...
        + x(nodes+2*edges+T).*basis_rt1_r(r,z,7)...
        + x(nodes+2*edges+triangles+T).*basis_rt1_r(r,z,8);
    
    approx_th = @(r,z) x(nodes+t_e(T,1)).*basis_rt1_th(r,z,1) + x(nodes+edges+t_e(T,1)).*basis_rt1_th(r,z,2)...
        + x(nodes+t_e(T,2)).*basis_rt1_th(r,z,3) + x(nodes+edges+t_e(T,2)).*basis_rt1_th(r,z,4)...
        + x(nodes+t_e(T,3)).*basis_rt1_th(r,z,5) + x(nodes+edges+t_e(T,3)).*basis_rt1_th(r,z,6)...
        + x(nodes+2*edges+T).*basis_rt1_th(r,z,7)...
        + x(nodes+2*edges+triangles+T).*basis_rt1_th(r,z,8)...
        + x(t(1,T)).*basis_p1_th(r,z,1)...
        + x(t(2,T)).*basis_p1_th(r,z,2)...
        + x(t(3,T)).*basis_p1_th(r,z,3); %2*edges+2*triangles+
    
    approx_z = @(r,z) x(nodes+t_e(T,1)).*basis_rt1_z(r,z,1)+x(nodes+edges+t_e(T,1)).*basis_rt1_z(r,z,2)...
        + x(nodes+t_e(T,2)).*basis_rt1_z(r,z,3) + x(nodes+edges+t_e(T,2)).*basis_rt1_z(r,z,4)...
        + x(nodes+t_e(T,3)).*basis_rt1_z(r,z,5) + x(nodes+edges+t_e(T,3)).*basis_rt1_z(r,z,6)...
        + x(nodes+2*edges+T).*basis_rt1_z(r,z,7)...
        + x(nodes+2*edges+triangles+T).*basis_rt1_z(r,z,8);
    
    % find L2 Error
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (u_vec_th(r,z) - approx_th(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z(r,z)).^2).*r;
 
    integral = integral + Wx'*feval(integrand,X,Y)*Wy;
end

err = sqrt(integral);

end

