function [err,grad_err,max_err] = errors_exact_weighted_HL_k_0_p2(p,t,p2,t2,basis,u_h,n,u,grad_u_r,grad_u_z)                                                             
% ERRORS_EXACT_WEIGHTED_HL_L_0_p2 - Calculate the errors of our solution u_h
% compared to the exact solution u.
%   Hodge Laplacian k = 0 case, P2
%
% Syntax:
%     [err,grad_err,max_err] = errors_exact_weighted_HL_k_0_p2(p,e,t,u_h_km1,u_h_k)
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
%         the pieceiwise basis function for the ith node in triangle T.
%     u_h - approximtate solution for u
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%     u - exact solution
%     grad_u_r - gradient(u) w.r.t r
%     grad_u_z - gradient(u) w.r.t z
%
% Outputs:
%    err - L2 error
%    grad_err - L2 gradient error
%    max_err - max error
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[~,mid_nodes] = size(p2);

integral = 0;
grad_integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
    
    [R,Z,Wr,Wz] = triquad(7, coordinates);

    b1 = basis(:,1,T);
    b2 = basis(:,2,T);
    b3 = basis(:,3,T);
    b4 = basis(:,4,T);
    b5 = basis(:,5,T);
    b6 = basis(:,6,T);
    
    n1 = t(1,T);
    n2 = t(2,T);
    n3 = t(3,T);
    n4 = t2(1,T);
    n5 = t2(2,T);
    n6 = t2(3,T);
    
    if n == 0
        % find L2 Error for u_h_km1 & u_h_k
        approx =@(r,z) u_h(n1).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h(n2).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h(n3).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h(n4 + nodes).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h(n5 + nodes).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h(n6 + nodes).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6));
    else
        % find L2 Error for u_h_km1 & u_h_k
        approx =@(r,z) u_h(n1).*(r./n).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h(n2).*(r./n).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h(n3).*(r./n).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h(n4 + nodes).*(r./n).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h(n5 + nodes).*(r./n).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h(n6 + nodes).*(r./n).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6));
    end
    
    integrand =@(r,z) ((u(r,z) - approx(r,z)).^2).*r;
    
    integral = integral + Wr'*feval(integrand,R,Z)*Wz;
    
    % find Grad L2 Error for u_h_km1 & u_h_k
    grad_approx_r =@(r,z) u_h(n1).*(2.*b1(1).*r + b1(2).*z ...
        + b1(4)) + u_h(n2).*(2.*b2(1).*r + b2(2).*z + b2(4)) ...
        + u_h(n3).*(2.*b3(1).*r + b3(2).*z + b3(4)) ...
        + u_h(n4 + nodes).*(2.*b4(1).*r + b4(2).*z + b4(4)) ...
        + u_h(n5 + nodes).*(2.*b5(1).*r + b5(2).*z + b5(4)) ...
        + u_h(n6 + nodes).*(2.*b6(1).*r + b6(2).*z + b6(4));
    grad_approx_z =@(r,z) u_h(n1).*(b1(2).*r + 2.*b1(3).*z ...
        + b1(5)) + u_h(n2).*(b2(2).*r + 2.*b2(3).*z + b2(5)) ...
        + u_h(n3).*(b3(2).*r + 2.*b3(3).*z + b3(5)) ...
        + u_h(n4 + nodes).*(b4(2).*r + 2.*b4(3).*z + b4(5)) ...
        + u_h(n5 + nodes).*(b5(2).*r + 2.*b5(3).*z + b5(5)) ...
        + u_h(n6 + nodes).*(b6(2).*r + 2.*b6(3).*z + b6(5));
        
    grad_integrand =@(r,z) ((grad_u_r(r,z) - grad_approx_r(r,z)).^2 ...
        + (grad_u_z(r,z) - grad_approx_z(r,z)).^2).*r;
    
    grad_integral = grad_integral + Wr'*feval(grad_integrand,R,Z)*Wz;
end

% find Max Error
max_err = -Inf;
for i = 1:nodes
    r = p(1,i);
    z = p(2,i);
    
    diff = abs(u(r,z) - u_h(i));
    
    if diff > max_err
        max_err = diff;
    end
end

for i = 1:mid_nodes
    r = p2(1,i);
    z = p2(2,i);
    
    diff = abs(u(r,z) - u_h(i + nodes));
    
    if diff > max_err
        max_err = diff;
    end
end

err = sqrt(integral);
grad_err = sqrt(grad_integral);

% end