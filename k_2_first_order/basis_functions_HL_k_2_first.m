function [basis_p1,basis_rt1,basis_p2,basis_nd1] = basis_functions_HL_k_2_first(p,t,p2,t2,ele,ed,new_ele)
% BASIS_FUNCTIONS_HL_K_2_first - Create a piecewise 
%   basis function for each node of a triangulation
%
% Syntax:
%     [basis_p1,basis_rt1,basis_p2,basis_nd1] = 
%           basis_functions_HL_k_2_first(p,t,p2,t2,ele,ed,new_ele)
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
%         each edge and triangle in each triangle. basis_rt1(:,i,T)
%         represents the pieceiwise basis function for the ith edge in
%         triangle T.
%     basis_p2 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_p2(:,i,T)
%         represents the pieceiwise basis function for the ith node or 
%         midpoint in triangle T.
%     basis_nd1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_nd1(:,i,T)
%         represents the ith pieceiwise basis function in triangle T.
%
% Dependencies:
%     basis_functions_weighted_HL_k_0_p1.m
%     basis_functions_rt1.m
%     basis_functions_weighted_HL_k_0_p2.m
%     basis_functions_nd1.m
%
% Author: Nicole Stock
% Date: Spring 2021
  
% get Ch1 Basis Functions basis functions
%   P1
basis_p1 = basis_functions_weighted_HL_k_0_p1(p,t);
%   RT1
basis_rt1 = basis_functions_rt1(p,t,ele,ed,new_ele);

% get Bh1 Basis Functions
%   P2
basis_p2 = basis_functions_weighted_HL_k_0_p2(p,t,p2,t2);
%   ND1
basis_nd1 = basis_functions_nd1(p,t,ele,ed,new_ele);

end