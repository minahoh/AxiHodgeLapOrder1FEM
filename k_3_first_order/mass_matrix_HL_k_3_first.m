function [mass_matrix] = mass_matrix_HL_k_3_first(p,t,ed,t_e,basis_p1,basis_rt1,n)
% MASS_MATRIX_HL_K_3_FIRST - Create mass matrix
%   Hodge Laplacian k = 3 case, First Order
%   (psi_i, psi_j)_r where {psi_k}k=1->(2Ne+2Nt+N) is the basis for Ch1
%   (Ch1 is the fourier first order Raviart Thomas space)
%
% Syntax:
%     A = mass_matrix_HL_k_3_first(p,t,ed,t_ed,basis_edges,basis_triangles,n)
%
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
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     mass_matrix - mass matrix
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*121);
j_vec = zeros(1,triangles*121);
s_vec = zeros(1,triangles*121);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
    
    % integrate for each pair of 2*edges + 2*triangles + nodes in the triangle
    for i = 1:11
        for j = i:11
            % get global i identity
            if i == 1
                global_i = nodes + t_e(T,1);
            elseif i == 2
                global_i = nodes + edges + t_e(T,1); 
            elseif i == 3
                global_i = nodes + t_e(T,2);
            elseif i == 4
                global_i = nodes + edges + t_e(T,2);
            elseif i == 5
                global_i = nodes + t_e(T,3);
            elseif i == 6
                global_i = nodes + edges + t_e(T,3);
            elseif i == 7
                global_i = nodes + 2*edges + T;
            elseif i == 8
                global_i = nodes + 2*edges + triangles + T; 
            else
                global_i = t(i-8,T);
            end
            % get global j identity
            if j == 1
                global_j = nodes + t_e(T,1);
            elseif j == 2
                global_j = nodes + edges + t_e(T,1); 
            elseif j == 3
                global_j = nodes + t_e(T,2);
            elseif j == 4
                global_j = nodes + edges + t_e(T,2);
            elseif j == 5
                global_j = nodes + t_e(T,3);
            elseif j == 6
                global_j = nodes + edges + t_e(T,3);
            elseif j == 7
                global_j = nodes + 2*edges + T;
            elseif j == 8
                global_j = nodes + 2*edges + triangles + T; 
            else
                global_j = t(j-8,T);
            end
            
            if i <= 8
                % Basis functions for edges and triangles using RT1
                phi_i_r = @(r,z) basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
                    + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z;
                phi_i_th = @(r,z) (1./n).*(basis_rt1(1,i,T).*r + basis_rt1(2,i,T).*z ...
                    + basis_rt1(3,i,T) + basis_rt1(7,i,T).*r.*r + basis_rt1(8,i,T).*r.*z);
                phi_i_z = @(r,z) basis_rt1(4,i,T).*r + basis_rt1(5,i,T).*z ...
                    + basis_rt1(6,i,T) + basis_rt1(7,i,T).*r.*z + basis_rt1(8,i,T).*z.*z;
            else
                % Basis functions for nodes using P1
                I = basis_p1(:,i-8,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                
                phi_i_r = @(r,z) 0;
                phi_i_th = @(r,z) (1./n).*(ai.*r.*r + bi.*r.*z + ci.*r);
                phi_i_z = @(r,z) 0;
            end
            if j <= 8
                % Basis functions for edges and triangles using RT1
                phi_j_r = @(r,z) basis_rt1(1,j,T).*r + basis_rt1(2,j,T).*z ...
                    + basis_rt1(3,j,T) + basis_rt1(7,j,T).*r.*r + basis_rt1(8,j,T).*r.*z;
                phi_j_th = @(r,z) (1./n).*(basis_rt1(1,j,T).*r + basis_rt1(2,j,T).*z ...
                    + basis_rt1(3,j,T) + basis_rt1(7,j,T).*r.*r + basis_rt1(8,j,T).*r.*z);
                phi_j_z = @(r,z) basis_rt1(4,j,T).*r + basis_rt1(5,j,T).*z ...
                    + basis_rt1(6,j,T) + basis_rt1(7,j,T).*r.*z + basis_rt1(8,j,T).*z.*z;
            else
                % Basis functions for nodes using P1
                J = basis_p1(:,j-8,T);
                aj = J(1);
                bj = J(2);
                cj = J(3);
                
                phi_j_r = @(r,z) 0;
                phi_j_th = @(r,z) (1./n).*(aj.*r.*r + bj.*r.*z + cj.*r);
                phi_j_z = @(r,z) 0;
            end

            % We assume n > 0 for all of these problems!

            integrand =@(r,z) (phi_i_r(r,z).*phi_j_r(r,z) ...
                + phi_i_th(r,z).*phi_j_th(r,z) ...
                + phi_i_z(r,z).*phi_j_z(r,z)).*r;
            
            Q = Wr'*feval(integrand,R,Z)*Wz;               
          
            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
            if global_j ~= global_i
                i_vec(index) = global_j;
                j_vec(index) = global_i;
                s_vec(index) = Q;
                index = index + 1;
            end
        end
    end
end
N = nodes + 2*edges + 2*triangles;
mass_matrix = sparse(i_vec,j_vec,s_vec,N,N);

end

