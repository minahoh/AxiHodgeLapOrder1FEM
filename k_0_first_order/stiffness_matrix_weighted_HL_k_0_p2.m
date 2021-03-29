function stiffness_matrix = stiffness_matrix_weighted_HL_k_0_p2(p,t,p2,t2,basis,n)
% STIFFNESS_MATRIX_WEIGHTED_HL_K_0_P2 - Create stiffness matrix with weight n
%
% Syntax:
%     A = stiffness_matrix_weighted_HL_k_0_p2(p,t,p2,t2,basis,k)
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
%     basis - a 6x6xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     stiffness_matrix - stiffness matrix
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*36);
j_vec = zeros(1,triangles*36);
s_vec = zeros(1,triangles*36);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates
        coordinates(N,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
    
    % integrate for each pair of nodes in the triangle
    for i = 1:6
        for j = i:6
            I = basis(:,i,T);
            J = basis(:,j,T);
            
            ai = I(1);
            bi = I(2);
            ci = I(3);
            di = I(4);
            ei = I(5);
            fi = I(6);
                        
            aj = J(1);
            bj = J(2);
            cj = J(3);
            dj = J(4);
            ej = J(5);
            fj = J(6);
            
            % integrate grad(I) * grad(J)
            % grad^k_rz(v) = [ partial_deriv_r(v)
            %                 (-k/r)*v
            %                 partial_deriv_z(v) ]
            % phi_k_i = (r/k)(ar^2 + brz + cz^2 + dr + ez + f)
            
            grad_i_r =@(r,z) (1./n).*(3.*ai.*r.^2 + 2.*bi.*r.*z + ci.*z.^2 ...
            + 2.*di.*r + ei.*z + fi);
            n_r_I =@(r,z) ai.*r.^2 + bi.*r.*z + ci.*z.^2 + di.*r ...
                + ei.*z + fi;
            grad_i_z =@(r,z) (1./n).*(bi.*r.^2 + 2.*ci.*r.*z + ei.*r);
            
            grad_j_r =@(r,z) (1./n).*(3.*aj.*r.^2 + 2.*bj.*r.*z + cj.*z.^2 ...
            + 2.*dj.*r + ej.*z + fj);
            n_r_J =@(r,z) aj.*r.^2 + bj.*r.*z + cj.*z.^2 + dj.*r ...
                + ej.*z + fj;
            grad_j_z =@(r,z) (1./n).*(bj.*r.^2 + 2.*cj.*r.*z + ej.*r);
            

            grad_integrand =@(r,z) (grad_i_r(r,z).*grad_j_r(r,z) ...
                + n_r_I(r,z).*n_r_J(r,z) ...
                + grad_i_z(r,z).*grad_j_z(r,z)).*r;
            Q = Wr'*feval(grad_integrand,R,Z)*Wz;               

            % if i or j == 4,5,6, they are midpoint nodes, so their global
            %   id is found in t2 added to the number of original nodes
            if i >= 4
                global_i = t2(i-3,T) + nodes;
            else
                global_i = t(i,T);
            end
            
            if j >=4
                global_j = t2(j-3,T) + nodes;
            else
                global_j = t(j,T);
            end
                        
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

[~,n2] = size(p2);
N = nodes + n2;
stiffness_matrix = sparse(i_vec,j_vec,s_vec,N,N);

% end