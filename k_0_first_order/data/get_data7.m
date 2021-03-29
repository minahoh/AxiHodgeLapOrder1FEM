function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data7()
%GET_U_DATA7 Get u and f function handles
%
%   u = r^(1/3)
% Author: Nicole Stock
% Date: Fall 2020

u = @(r,z) r.^(1./3);
grad_u_r =@(r,z) (1./3).*r.^(-2./3);
grad_u_z =@(r,z) 0;

f = u;
grad_f_r = grad_u_r;
grad_f_z = grad_u_z;
end

