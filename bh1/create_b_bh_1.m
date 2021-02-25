function b = create_b_bh_1(p,t,p2,t2,ed,t_e,basis_p2,basis_nd1,f_vec_r,f_vec_th,f_vec_z,n)
% CREATE_B_BH_1 - Create vector b such that 
%     stiffness_matrix * solution = b.
%
% Syntax:
%     b = create_b_bh_1(p,t,ed,t_ed,basis,f_vec_r,f_vec_z)
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
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,T) represents 
%         the pieceiwise basis function for the ith node in triangle T. 
%     f_vec_r - given vector r component
%     f_vec_z - given vector z component
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,nodes] = size(p);
[~,midpoints] = size(p2);
[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*14);
j_vec = ones(1,triangles*14);
s_vec = zeros(1,triangles*14);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    % for each edge in the triangle
    for i = 1:14
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
        
        integrand =@(r,z) (f_vec_r(r,z).*phi_i_r(r,z) ...
            + f_vec_th(r,z).*phi_i_th(r,z) ...
            + f_vec_z(r,z).*phi_i_z(r,z)).*r;
        
        b_i = Wx'*feval(integrand,X,Y)*Wy;
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end
N = nodes + midpoints + 2*edges + 2*triangles;
b = sparse(i_vec,j_vec,s_vec,N,1);