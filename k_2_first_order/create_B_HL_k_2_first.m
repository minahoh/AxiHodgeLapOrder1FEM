function B = create_B_HL_k_2_first(p,t,p2,t2,ed,t_e,basis_p1,basis_rt1,basis_p2,basis_nd1,n)
% CREATE_B_HL_K_2_first - Create B matrix
%   Hodge Laplacian k = 2 case, first order
%   (curl_rz^n(zeta_i), psi_i)_r where {zeta_i} is the basis for Bh1 and
%     {psi_j}t is the basis for Ch1
%   (Bh1 is the weighted fourier ND1 and P2 space)
%   (Ch1 is the weighted fourier RT1 and P1 space)
%
% Syntax:
%     B = create_B_HL_k_2_first(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,n)
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
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     B - B matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,midpoints] = size(p2);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*154);
j_vec = zeros(1,triangles*154);
s_vec = zeros(1,triangles*154);
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
    
    % integrate for each pair of (nodees,midpoints,2*edges,2*triangles) 
    %   x (2*edges,2*triangles,nodes) in the triangle
    for i = 1:14
        for j = 1:11
            % get global i identity
            if i <= 3
                % node
                global_i = 2*edges + 2*triangles + t(i,T);
            elseif i <= 6
                % midpoint
                global_i = 2*edges + 2*triangles + nodes + t2(i-3,T);
            elseif i == 7
                global_i = t_e(T,1);
            elseif i == 8
                global_i = edges + t_e(T,1); 
            elseif i == 9
                global_i = t_e(T,2);
            elseif i == 10
                global_i = edges + t_e(T,2);
            elseif i == 11
                global_i = t_e(T,3);
            elseif i == 12
                global_i = edges + t_e(T,3);
            elseif i == 13
                global_i = 2*edges + T;
            else
                % i == 14
                global_i = 2*edges + triangles + T; 
            end
            % get global j identity
            if j == 1
                global_j = t_e(T,1);
            elseif j == 2
                global_j = edges + t_e(T,1); 
            elseif j == 3
                global_j = t_e(T,2);
            elseif j == 4
                global_j = edges + t_e(T,2);
            elseif j == 5
                global_j = t_e(T,3);
            elseif j == 6
                global_j = edges + t_e(T,3);
            elseif j == 7
                global_j = 2*edges + T;
            elseif j == 8
                global_j = 2*edges + triangles + T; 
            else
                global_j = 2*edges + 2*triangles + t(j-8,T);
            end
            
            if i <=6
                % P2 Basis Function
                I = basis_p2(:,i,T);

                ai = I(1);
                bi = I(2);
                ci = I(3);
                di = I(4);
                ei = I(5);
                %fi = I(6);

                curl_i_r = @(r,z) -(bi.*r + 2.*ci.*z + ei);
                curl_i_th = @(r,z) (-1./n).*(bi.*r + 2.*ci.*z + ei);
                curl_i_z = @(r,z) 2.*ai.*r + bi.*z + di;
            else
                % ND1 Basis Function
                i2 = i - 6;
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
            
            integrand =@(r,z) (curl_i_r(r,z).*phi_j_r(r,z) ...
                + curl_i_th(r,z).*phi_j_th(r,z) ...
                + curl_i_z(r,z).*phi_j_z(r,z)).*r;

            Q = Wr'*feval(integrand,R,Z)*Wz;               

            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
        end
    end
end

N1 = nodes + midpoints + 2*edges + 2*triangles;
N2 = nodes + 2*edges + 2*triangles;
B = sparse(i_vec,j_vec,s_vec,N1,N2);
% end