function [err,grad_err] = errors_no_exact_HL_k_3_first(p,t,basis_p1,u_h,u_h_km1)
% ERRORS_NO_EXACT_HL_K_3_FIRST - Calculate the errors of our solution x
% compared to the exact solution u.
%   Hodge Laplacian k = 3 case, first order
%
% Syntax:
%     [err_z,err_p] = 
%         errors_no_exact_HL_k_3_first(p,t,basis_p1,u_h,u_h_km1)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     basis_p1 - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_vertices(:,i,T)
%         represents the pieceiwise basis function for the ith edge in 
%         triangle T.
%     u_h - approximated solution of p
%     u_h_km1 - approximate solution for mesh level k-1
%
% Outputs:
%    err - L2 error for u approximation
%
% Author: Nicole Stock
% Date: Spring 2021

[~,triangles] = size(t);

integral = 0;
grad_integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates
        coordinates(N,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    % find L2 Error for p
    basis_p1_fn =@(r,z,i) basis_p1(1,i,T).*r ... 
        + basis_p1(2,i,T).*z + basis_p1(3,i,T);
    
    approx = @(r,z) u_h(t(1,T)).*basis_p1_fn(r,z,1)...
        + u_h(t(2,T)).*basis_p1_fn(r,z,2)...
        + u_h(t(3,T)).*basis_p1_fn(r,z,3);
    approx_km1 = @(r,z) u_h_km1(t(1,T)).*basis_p1_fn(r,z,1)...
        + u_h_km1(t(2,T)).*basis_p1_fn(r,z,2)...
        + u_h_km1(t(3,T)).*basis_p1_fn(r,z,3);
    
    integrand =@(r,z) ((approx_km1(r,z) - approx(r,z)).^2).*r;
    
    integral = integral + Wx'*feval(integrand,X,Y)*Wy; 
    
    % find Grad L2 Error
    b1 = basis_p1(:,1,T);
    b2 = basis_p1(:,2,T);
    b3 = basis_p1(:,3,T);
    
    n1 = t(1,T);
    n2 = t(2,T);
    n3 = t(3,T);
    
    grad_approx_x_km1 = u_h_km1(n1).*b1(1) + u_h_km1(n2).*b2(1) ...
        + u_h_km1(n3).*b3(1);
    grad_approx_y_km1 = u_h_km1(n1).*b1(2) + u_h_km1(n2).*b2(2) ...
        + u_h_km1(n3).*b3(2);
    
    grad_approx_x_k = u_h(n1).*b1(1) + u_h(n2).*b2(1) ...
        + u_h(n3).*b3(1);
    grad_approx_y_k = u_h(n1).*b1(2) + u_h(n2).*b2(2) ...
        + u_h(n3).*b3(2);  

    grad_integrand =@(x,y) ((grad_approx_x_km1 - grad_approx_x_k).^2 ...
        + (grad_approx_y_km1 - grad_approx_y_k).^2).*x;

    grad_integral = grad_integral + Wx'*feval(grad_integrand,X,Y)*Wy;
end

err = sqrt(integral);
grad_err = sqrt(grad_integral);

% end