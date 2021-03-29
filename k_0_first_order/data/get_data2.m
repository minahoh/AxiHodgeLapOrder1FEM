function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data2()
%GET_U_DATA1 Get f, grad_rz^1(f) wrt r grad_rz^1(f) wrt z, u, grad_rz^1(u) 
%   wrt r and grad_rz^1(u) wrt z
%   data2
%   u = r^3 - (3/2)r^2
%   f = -8r + 9/2
% Author: Nicole Stock
% Date: Fall 2020

f = @(r,z) -8.*r + (9./2);
grad_f_r =@(r,z) -8;
grad_f_z =@(r,z) 0;

u = @(r,z) r.^3 - (3./2).*r.^2;
grad_u_r =@(r,z) 3.*r.^2 - 3.*r;
grad_u_z =@(r,z) 0;
end