function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data6()
%GET_U_DATA6 Get u and f function handles
%
%   u = r*sin(z)
% Author: Nicole Stock
% Date: Fall 2020

u = @(r,z) r.*sin(z);
grad_u_r =@(r,z) sin(z);
grad_u_z =@(r,z) r.*cos(z);

f = u;
grad_f_r = grad_u_r;
grad_f_z = grad_u_z;
end

