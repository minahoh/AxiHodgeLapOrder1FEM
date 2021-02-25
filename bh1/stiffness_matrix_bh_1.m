function stiffness_matrix = stiffness_matrix_bh_1(p,t,p2,t2,ed,t_e,basis_p2,basis_nd1,n)
% STIFFNESS_MATRIX_BH_1 - Create stiffness matrix
%
% Syntax:
%     A = stiffness_matrix_bh_1(p,t,ed,t_ed,basis)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_e - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
%     basis_p2 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_p2(:,i,T)
%         represents the pieceiwise basis function for the ith node or 
%         midpoint in triangle T.
%     basis_nd1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_nd1(:,i,T)
%         represents the ith pieceiwise basis function in triangle T.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     stiffness_matrix - stiffness matrix
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,midpoints] = size(p2);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*196);
j_vec = zeros(1,triangles*196);
s_vec = zeros(1,triangles*196);
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
    
    % integrate for each pair of edges in the triangle
    for i = 1:14
        for j = i:14
            % get global i identity
            if i <= 3
                % node
                global_i = t(i,T);
            elseif i <= 6
                % midpiint
                global_i = nodes + t2(i-3,T);
            elseif i == 7
                global_i = nodes + midpoints + t_e(T,1);
            elseif i == 8
                global_i = nodes + midpoints + edges + t_e(T,1); 
            elseif i == 9
                global_i = nodes + midpoints + t_e(T,2);
            elseif i == 10
                global_i = nodes + midpoints + edges + t_e(T,2);
            elseif i == 11
                global_i = nodes + midpoints + t_e(T,3);
            elseif i == 12
                global_i = nodes + midpoints + edges + t_e(T,3);
            elseif i == 13
                global_i = nodes + midpoints + 2*edges + T;
            else
                % i == 14
                global_i = nodes + midpoints + 2*edges + triangles + T; 
            end
            % get global j identity
            if j <= 3
                % node
                global_j = t(j,T);
            elseif j <= 6
                % midpiint
                global_j = nodes + t2(j-3,T);
            elseif j == 7
                global_j = nodes + midpoints + t_e(T,1);
            elseif j == 8
                global_j = nodes + midpoints + edges + t_e(T,1); 
            elseif j == 9
                global_j = nodes + midpoints + t_e(T,2);
            elseif j == 10
                global_j = nodes + midpoints + edges + t_e(T,2);
            elseif j == 11
                global_j = nodes + midpoints + t_e(T,3);
            elseif j == 12
                global_j = nodes + midpoints + edges + t_e(T,3);
            elseif j == 13
                global_j = nodes + midpoints + 2*edges + T;
            else
                % j == 14
                global_j = nodes + midpoints + 2*edges + triangles + T; 
            end
            
            if i <=6
                % P2 Basis Function
                I = basis_p2(:,i,T);

                ai = I(1);
                bi = I(2);
                ci = I(3);
                di = I(4);
                ei = I(5);
                fi = I(6);
                
                phi_i_r = @(r,z) (-1./n).*(ai.*r.^2 + bi.*r.*z ...
                    + ci.*z.^2 + di.*r + ei.*z + fi);
                phi_i_th = @(r,z) ai.*r.^2 + bi.*r.*z + ci.*z.^2 ...
                    + di.*r + ei.*z + fi;
                phi_i_z = @(r,z) 0;
                
                curl_i_r = @(r,z) -(bi.*r + 2.*ci.*z + ei);
                curl_i_th = @(r,z) (-1./n).*(bi.*r + 2.*ci.*z + ei);
                curl_i_z = @(r,z) 2.*ai.*r + bi.*z + di;
            else
                % ND1 Basis Function
                i2 = i - 6;
                phi_i_r = @(r,z) (1./n).*(basis_nd1(1,i2,T).*r.^2 ...
                    + basis_nd1(2,i2,T).*r.*z + basis_nd1(3,i2,T).*r ...
                    - basis_nd1(7,i2,T).*r.^2.*z - basis_nd1(8,i2,T).*r.*z.^2);
                phi_i_th = @(r,z) 0;
                phi_i_z = @(r,z) (1./n).*(basis_nd1(4,i2,T).*r.^2 ...
                    + basis_nd1(5,i2,T).*r.*z + basis_nd1(6,i2,T).*r ...
                    + basis_nd1(7,i2,T).*r.^3 + basis_nd1(8,i2,T).*r.^2.*z);
                
                curl_i_r = @(r,z) -(basis_nd1(4,i2,T).*r ...
                    + basis_nd1(5,i2,T).*z + basis_nd1(6,i2,T) ...
                    + basis_nd1(7,i2,T).*r.^2 + basis_nd1(8,i2,T).*r.*z);
                curl_i_th = @(r,z) (-1./n).*(-basis_nd1(2,i2,T).*r ...
                    + 2.*basis_nd1(4,i2,T).*r + basis_nd1(5,i2,T).*z ...
                    + basis_nd1(6,i2,T) + 4.*basis_nd1(7,i2,T).*r.^2 ...
                    + 4.*basis_nd1(8,i2,T).*r.*z);
                curl_i_z = @(r,z) basis_nd1(1,i2,T).*r ...
                    + basis_nd1(2,i2,T).*z + basis_nd1(3,i2,T) ...
                    - basis_nd1(7,i2,T).*r.*z - basis_nd1(8,i2,T).*z.^2;  
            end
            
            if j <=6
                % P2 Basis Function
                J = basis_p2(:,j,T);

                aj = J(1);
                bj = J(2);
                cj = J(3);
                dj = J(4);
                ej = J(5);
                fj = J(6);
                
                phi_j_r = @(r,z) (-1./n).*(aj.*r.^2 + bj.*r.*z ...
                    + cj.*z.^2 + dj.*r + ej.*z + fj);
                phi_j_th = @(r,z) aj.*r.^2 + bj.*r.*z + cj.*z.^2 ...
                    + dj.*r + ej.*z + fj;
                phi_j_z = @(r,z) 0;
                
                curl_j_r = @(r,z) -(bj.*r + 2.*cj.*z + ej);
                curl_j_th = @(r,z) (-1./n).*(bj.*r + 2.*cj.*z + ej);
                curl_j_z = @(r,z) 2.*aj.*r + bj.*z + dj;
            else
                % ND1 Basis Function
                j2 = j - 6;
                phi_j_r = @(r,z) (1./n).*(basis_nd1(1,j2,T).*r.^2 ...
                    + basis_nd1(2,j2,T).*r.*z + basis_nd1(3,j2,T).*r ...
                    - basis_nd1(7,j2,T).*r.^2.*z - basis_nd1(8,j2,T).*r.*z.^2);
                phi_j_th = @(r,z) 0;
                phi_j_z = @(r,z) (1./n).*(basis_nd1(4,j2,T).*r.^2 ...
                    + basis_nd1(5,j2,T).*r.*z + basis_nd1(6,j2,T).*r ...
                    + basis_nd1(7,j2,T).*r.^3 + basis_nd1(8,j2,T).*r.^2.*z);
                
                curl_j_r = @(r,z) -(basis_nd1(4,j2,T).*r ...
                    + basis_nd1(5,j2,T).*z + basis_nd1(6,j2,T) ...
                    + basis_nd1(7,j2,T).*r.^2 + basis_nd1(8,j2,T).*r.*z);
                curl_j_th = @(r,z) (-1./n).*(-basis_nd1(2,j2,T).*r ...
                    + 2.*basis_nd1(4,j2,T).*r + basis_nd1(5,j2,T).*z ...
                    + basis_nd1(6,j2,T) + 4.*basis_nd1(7,j2,T).*r.^2 ...
                    + 4.*basis_nd1(8,j2,T).*r.*z);
                curl_j_z = @(r,z) basis_nd1(1,j2,T).*r ...
                    + basis_nd1(2,j2,T).*z + basis_nd1(3,j2,T) ...
                    - basis_nd1(7,j2,T).*r.*z - basis_nd1(8,j2,T).*z.^2;
            end
                        
            integrand =@(r,z) (phi_i_r(r,z).*phi_j_r(r,z) ...
                + phi_i_th(r,z).*phi_j_th(r,z) ...
                + phi_i_z(r,z).*phi_j_z(r,z) ...
                + curl_i_r(r,z).*curl_j_r(r,z) ...
                + curl_i_th(r,z).*curl_j_th(r,z) ...
                + curl_i_z(r,z).*curl_j_z(r,z)).*r;
            
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

N = nodes + midpoints + 2*edges + 2*triangles;
stiffness_matrix = sparse(i_vec,j_vec,s_vec,N,N);

% end