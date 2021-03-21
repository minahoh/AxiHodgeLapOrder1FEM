function mass_matrix = mass_matrix_HL_k_1_first(p,t,p2,t2,basis_p2,n)
% MASS_MATRIX_HL_K_1_FIRST - Create mass matrix
%   Hodge Laplacian k = 1 case, first order
%   (phi_i, phi_j)_r where {phi_k} is the basis for Ah1
%   (Ah1 is the weighted P2 space)
%
% Syntax:
%     M = mass_matrix_HL_k_1_first(p,t,ed,t_ed,basis_nodes,basis_edges,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     new_ele - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
%     basis_p2 - a 6x6xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     n - Fourier mdoe
%
% Outputs:
%     mass_matrix - mass matrix used to solve system of equations to
%         approximate solution
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,midpoints] = size(p2);
[~,triangles] = size(t);
i_vec = zeros(1,triangles*36);
j_vec = zeros(1,triangles*36);
s_vec = zeros(1,triangles*36);
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
    
    % integrate for each pair of vertices in the triangle
    for i = 1:6
        for j = i:6
            % get global i identity
            if i <= 3
                % node
                global_i = t(i,T);
            elseif i <= 6
                % midpoint
                global_i = nodes + t2(i-3,T);
            end
            % get global j identity
            if j <= 3
                % node
                global_j = t(j,T);
            else
                % midpiint
                global_j = nodes + t2(j-3,T);
            end
            
            I = basis_p2(:,i,T);

            ai = I(1);
            bi = I(2);
            ci = I(3);
            di = I(4);
            ei = I(5);
            fi = I(6);
            phi_i =@(r,z) (1./n).*(ai.*r.^3 + bi.*r.^2.*z + ci.*r.*z.^2 ...
                    + di.*r.^2 + ei.*r.*z + fi.*r);

            J = basis_p2(:,j,T);

            aj = J(1);
            bj = J(2);
            cj = J(3);
            dj = J(4);
            ej = J(5);
            fj = J(6);
            phi_j =@(r,z) (1./n).*(aj.*r.^3 + bj.*r.^2.*z + cj.*r.*z.^2 ...
                    + dj.*r.^2 + ej.*r.*z + fj.*r);
 
            integrand =@(r,z) phi_i(r,z).*phi_j(r,z).*r;
            
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

N = nodes+midpoints;
mass_matrix = sparse(i_vec,j_vec,s_vec,N,N);
% end