function [basis_p1,basis_rt1] = basis_functions_HL_k_3_first(p,t,ele,ed,new_ele)
%BASIS_FUNCTIONS_RT1 Create a piecewise basis function for 
%   each edge and triangle of a triangulation
%   Hodge Laplacian k = 3 case, Ch1 & Dh1, first order
%
% Syntax:
%     [basis_rt1,basis_p1] = basis_functions_HL_k_3_first(p,t,ele,ed,new_ele)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ele - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     new_ele - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
%
% Outputs:
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     basis_rt1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%
% Author: Nicole Stock
% Date: Spring 2021

% get Ch1 Basis Functions
%   RT1 basis functions
basis_rt1 = basis_functions_rt1(p,t,ele,ed,new_ele);
%   P1 basis functions
basis_p1 = basis_functions_weighted_HL_k_0_p1(p,t);

end

