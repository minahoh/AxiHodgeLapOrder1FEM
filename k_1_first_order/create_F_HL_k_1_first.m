function F = create_F_HL_k_1_first(p,t,p2,t2,ed,t_e,basis_p2,basis_nd1,f_r,f_th,f_z,n)
% CREATE_F_HL_K_1_FIRST - Create F matrix
%   Hodge Laplacian k = 1 case, first order
%   (f, zeta_i)_r where {zeta_j} is the basis for Bh1
%   (Bh1 is the weighted fourier ND1 and P2 space)
%
% Syntax:
%     F = create_F_HL_k_1_first(p,t,p2,t2,ed,t_e,basis_p2,basis_nd1,f_r,f_th,f_z,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     f_vec_r - given function r component
%     f_vec_th - given function theta component
%     f_vec_z - given function z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     F - F matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,midpoints] = size(p2);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*14);
j_vec = ones(1,triangles*14);
s_vec = zeros(1,triangles*14);
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
    
    % integrate for each pair of basis functions
    for i = 1:14
        % get global i identity
        if i <= 3
            % node
            global_i= t(i,T);
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
        end

        integrand =@(r,z) (phi_i_r(r,z).*f_r(r,z) ...
            + phi_i_th(r,z).*f_th(r,z) ...
            + phi_i_z(r,z).*f_z(r,z)).*r;

        Q = Wr'*feval(integrand,R,Z)*Wz;               

        i_vec(index) = global_i;
        s_vec(index) = Q;
        index = index + 1;
    end
end

N = nodes + midpoints + 2*edges + 2*triangles;
F = sparse(i_vec,j_vec,s_vec,N,1);

% end