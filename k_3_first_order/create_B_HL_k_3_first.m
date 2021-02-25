function B = create_B_HL_k_3_first(p,t,ed,t_e,basis_p1,basis_rt1)
% CREATE_B_HL_K_3_FIRST - Create mass matrix
%   Hodge Laplacian k = 3 case, first order
%   (chi_i, div_rz^n(psi_j))_r where {psi_i}i=1->2Ne+2Nt+N is the basis for
%     Ch1 and {chi_j}j=1->N is the basis for Dh1
%   (Ch1 is the fourier first order Raviart Thomas space)
%   (Dh1 is the piecewise lineear space)
%
% Syntax:
%     B = create_B_HL_k_3_first(p,t,ed,t_e,basis)
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
%     basis_vertices - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%
% Outputs:
%     B - B matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*33);
j_vec = zeros(1,triangles*33);
s_vec = zeros(1,triangles*33);
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
        for j = 1:3
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
                global_i = t(i-8,T); % 2*edges + 2*triangles + 
            end
            
            if i <= 8
                % Basis functions for edges and triangles in Ch1
                div_i = @(r,z) basis_rt1(1,i,T) + basis_rt1(5,i,T) ...
                    + 3.*basis_rt1(7,i,T).*r + 3.*basis_rt1(8,i,T).*z;
            else
                % Basis functions for nodes in Ch1
                I = basis_p1(:,i-8,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                div_i = @(r,z) -(ai.*r + bi.*z + ci);
            end
            
            J = basis_p1(:,j,T);
            aj = J(1);
            bj = J(2);
            cj = J(3);
            
            phi_j = @(r,z) aj.*r + bj.*z + cj;
            
            global_j = t(j,T);

            integrand =@(r,z) (div_i(r,z).*phi_j(r,z)).*r;

            Q = Wr'*feval(integrand,R,Z)*Wz;               

            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
        end
    end
end

N1 = 2*edges + 2*triangles + nodes;
N2 = nodes;
B = sparse(i_vec,j_vec,s_vec,N1,N2);

end

