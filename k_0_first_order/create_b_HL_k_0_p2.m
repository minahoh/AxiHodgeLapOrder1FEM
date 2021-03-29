function b = create_b_HL_k_0_p2(p,t,p2,t2,basis,f_fn,grad_f_r,grad_f_z,n)
% CREATE_B_HL_K_0_P2 - Create vector b such that
%     stiffness_matrix * solution = b.
%   Hodge Laplacian k = 0 case, P2
%
% Syntax:
%     b = create_b_HL_p2(p,t,p2,t2,basis,f_fn,grad_f_r,grad_f_z,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     f - given function
%     grad_f_r - gradient(f) w.r.t r
%     grad_f_z - gradient(f) w.r.t z
%     n - Hodge Laplacian on Axisymmetrix Domain and its Discretization
%     weight
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*6);
j_vec = ones(1,triangles*6);
s_vec = zeros(1,triangles*6);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    for i = 1:6
        
        I = basis(:,i,T);
        a = I(1);
        b = I(2);
        c = I(3);
        d = I(4);
        e = I(5);
        f = I(6);
        
        % Calculate the gradient_rz^n* of basis function I
        grad_I_r =@(r,z) (1./n).*(3.*a.*r.^2 + 2.*b.*r.*z + c.*z.^2 ...
            + 2.*d.*r + e.*z + f);
        n_r_I =@(r,z) a.*r.^2 + b.*r.*z + c.*z.^2 + d.*r + e.*z + f;
        grad_I_z =@(r,z) (1./n).*(b.*r.^2 + 2.*c.*r.*z + e.*r);
        
        integrand =@(r,z) (grad_f_r(r,z).*grad_I_r(r,z) ...
            + (n./r).*f_fn(r,z).*n_r_I(r,z) ...
            + grad_f_z(r,z).*grad_I_z(r,z)).*r;
        
        % Integrate and solve
        b_i = Wx'*feval(integrand,X,Y)*Wy;
        
        if i >= 4
            global_i = t2(i-3,T) + nodes;
        else
            global_i = t(i,T);
        end
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

[~,n2] = size(p2);
b = sparse(i_vec,j_vec,s_vec,nodes+n2,1);