function basis = basis_functions_weighted_HL_k_0_p2(p,t,p2,t2)
% BASIS_FUNCTIONS_WEIGHTED_HL_K_0_P2 - Create a piecewise basis function for each
% node of a triangulation with weight k
%   Hodge Laplacian k = 0 case, P2
%
% Syntax:
%     basis = basis_functions_weighted_HL_k_0_p2(p,t,p2,t2)
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
%
% Outputs:
%     basis - a matrix representing piece-wise basis functions for each node
%         in each triangle. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith node in triangle T. 
%
% Author: Nicole Stock
% Date: Fall 2020
    
[~,triangles] = size(t);
basis = zeros(6,6,triangles);

for T = 1:triangles
    poly = zeros(6,6);
    
    % for each row (node)
    for i = 1:3
        % 'regular' node
        node = t(i, T);
        % get r,z coordinates
        r = p(1,node);
        z = p(2,node);
        
        % fill in polynomial
        poly(i,1) = r.^2;
        poly(i,2) = r.*z;
        poly(i,3) = z.^2;
        poly(i,4) = r;
        poly(i,5) = z;
        poly(i,6) = 1;
        
        % midpoint node
        node = t2(i,T);
        % get r,z coordinates
        r = p2(1,node);
        z = p2(2,node);
        % fill in polynomial
        poly(i+3,1) = r.^2;
        poly(i+3,2) = r.*z;
        poly(i+3,3) = z.^2;
        poly(i+3,4) = r;
        poly(i+3,5) = z;
        poly(i+3,6) = 1;
    end
    
    I = eye(6);
    basis(:,:,T) = poly\I;
end

% end
