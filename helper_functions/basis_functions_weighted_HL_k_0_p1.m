function basis = basis_functions_weighted_HL_k_0_p1(p,t)
% BASIS_FUNCTIONS_WEIGHTED_HL_K_0_P1 - Create a piecewise basis function for each
% node of a triangulation with weight k
%   Hodge Laplacian k = 0 case, P1
%
% Syntax:
%     basis = basis_functions_weighted_HL_k_0_p1(p,t,p2,t2)
% 
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%
% Outputs:
%     basis - a matrix representing piece-wise basis functions for each node
%         in each triangle. basis(:,i,T) represents the pieceiwise basis 
%         function for the ith node in triangle T. 
%
% Author: Nicole Stock
% Date: Spring 2020
    
[~,triangles] = size(t);
basis = zeros(3,3,triangles);

for T = 1:triangles
    poly = zeros(3,3);
    
    % for each row (node)
    for i = 1:3
        % 'regular' node
        node = t(i, T);
        % get r,z coordinates
        r = p(1,node);
        z = p(2,node);
        
        % fill in polynomial
        poly(i,1) = r;
        poly(i,2) = z;
        poly(i,3) = 1;
    end
    
    I = eye(3);
    basis(:,:,T) = poly\I;
end

% end
