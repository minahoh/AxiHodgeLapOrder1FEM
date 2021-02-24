function [b_vector] = create_b_ch_1(p,t,ed,t_e,basis_vertices,basis_rt1,f_vec_r,f_vec_th,f_vec_z,n)
% CREATE_B_WEIGHTED_CH_1 - Create vector b such that 
%     stiffness_matrix * solution = b.
%
% Syntax:
%     b = create_b_ch_1(p,t,ed,t_e,basis_vertices,basis_rt1,f_vec_r,f_vec_th,f_vec_z,n)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
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
%     f_vec_r - given vector r component
%     f_vec_th - given vector theta component
%     f_vec_z - given vector z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*11);
j_vec = ones(1,triangles*11);
s_vec = zeros(1,triangles*11);
index = 1;

for T = 1:triangles
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for s = 1:3
        node = t(s,T);
        % get x,y coordinates of triangle
        coordinates(s,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
    
    % integrate for each pair of edges + 1 in the triangle
    for i = 1:11
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
            I = basis_vertices(:,i-8,T);
            ai = I(1);
            bi = I(2);
            ci = I(3);

            phi_i_r = @(r,z) 0;
            phi_i_th = @(r,z) (1./n).*(ai.*r.*r + bi.*r.*z + ci.*r);
            phi_i_z = @(r,z) 0;
        end
        
        integrand =@(r,z) (f_vec_r(r,z).*phi_i_r(r,z) ...
            + f_vec_th(r,z).*phi_i_th(r,z) ...
            + f_vec_z(r,z).*phi_i_z(r,z)).*r;

        Q = Wr'*feval(integrand,R,Z)*Wz;               

        i_vec(index) = global_i;
        s_vec(index) = Q;
        index = index + 1;
    end
    
end

N = nodes + 2*edges + 2*triangles;
b_vector = sparse(i_vec,j_vec,s_vec,N,1);

end



