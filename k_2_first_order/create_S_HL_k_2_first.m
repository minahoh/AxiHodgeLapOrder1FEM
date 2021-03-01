function S = create_S_HL_k_2_first(p,t,ed,t_e,basis_p1,basis_rt1)
% CREATE_S_HL_K_2_FIRST - Create S matrix
%   Hodge Laplacian k = 2 case, first order
%   (div_rz^n(psi_i), div_rz^n(psi_j))_r where {psi_j} is the basis for Ch1
%   (Ch1 is the weighted fourier RT1 and P1 space)
%
% Syntax:
%     S = create_S_HL_k_2_first(p,t,ed,t_ed,basis_RT_edges,basis_triangles)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T)
%         represents the pieceiwise basis function for the ith edge in
%         triangle T.
%
% Outputs:
%     S - S matrix used to solve system of equations to approximate
%         solution
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
    
    % integrate for each pair of edges + 1 in the triangle
    for i = 1:11
        for j = i:11
            % get global i identity
            if i == 1
                global_i = t_e(T,1);
            elseif i == 2
                global_i = edges + t_e(T,1); 
            elseif i == 3
                global_i = t_e(T,2);
            elseif i == 4
                global_i = edges + t_e(T,2);
            elseif i == 5
                global_i = t_e(T,3);
            elseif i == 6
                global_i = edges + t_e(T,3);
            elseif i == 7
                global_i = 2*edges + T;
            elseif i == 8
                global_i = 2*edges + triangles + T; 
            else
                global_i = 2*edges + 2*triangles + t(i-8,T);
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

            if i <= 8
                % Basis functions for edges and triangles using RT1
                 div_i = @(r,z) basis_rt1(1,i,T) + basis_rt1(5,i,T) ...
                    + 3.*basis_rt1(7,i,T).*r + 3.*basis_rt1(8,i,T).*z;
            else
                % Basis functions for nodes using P1
                I = basis_p1(:,i-8,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                
                div_i = @(r,z) -(ai.*r + bi.*z + ci);
            end
            if j <= 8
                % Basis functions for edges and triangles using RT1
                div_j = @(r,z) basis_rt1(1,j,T) + basis_rt1(5,j,T) ...
                    + 3.*basis_rt1(7,j,T).*r + 3.*basis_rt1(8,j,T).*z;
            else
                % Basis functions for nodes using P1
                J = basis_p1(:,j-8,T);
                aj = J(1);
                bj = J(2);
                cj = J(3);
                
                div_j = @(r,z) -(aj.*r + bj.*z + cj);
            end
            
            integrand =@(r,z) (div_i(r,z).*div_j(r,z)).*r;

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
S = sparse(i_vec,j_vec,s_vec,N,N);

% end