function F = create_F_HL_k_3_first(p,t,basis_p1,f)
% CREATE_F_HL_K_3_FIRST - Create mass matrix
%   Hodge Laplacian k = 3 case, first order
%   (f, chi_i)_r where {chi_j}j=1->N is the basis for Dh1
%   (Dh1 is the piecewise linear space)
%
% Syntax:
%      F = create_F_HL_k_3_first(p,t,basis_vertices,f)
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
%     f - given function
%
% Outputs:
%     F - F matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Spring 2021

[~,nodes] = size(p);
[~,triangles] = size(t);
i_vec = zeros(1,triangles*3);
j_vec = ones(1,triangles*3);
s_vec = zeros(1,triangles*3);
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
    
    for i = 1:3
        I = basis_p1(:,i,T);
        ai = I(1);
        bi = I(2);
        ci = I(3);

        phi_i = @(r,z) ai.*r + bi.*z + ci;

        global_i = t(i,T);
        
        integrand =@(r,z) (f(r,z).*phi_i(r,z)).*r;

        Q = Wr'*feval(integrand,R,Z)*Wz;               

        i_vec(index) = global_i;
        s_vec(index) = Q;
        index = index + 1;
    end
end

F = sparse(i_vec,j_vec,s_vec,nodes,1);