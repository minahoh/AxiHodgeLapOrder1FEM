function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data1()
%GET_DATA1 Get f, grad_rz^1(f) wrt r, grad_rz^1(f) wrt z, u, grad_rz^1(u) 
%   wrt r and grad_rz^1(u) wrt z for data1
%   data1
%   u = (1/3)r^3 - (1/2)r^2
%   f = -(8/3)r + (5/2)
% Author: Nicole Stock
% Date: Fall 2020

f = @(r,z) -(8./3).*r + (3./2);
grad_f_r =@(r,z) -8./3;
grad_f_z =@(r,z) 0;

u = @(r,z) (1./3).*r.^3 - (1./2).*r.^2;
grad_u_r =@(r,z) r.^2 - r;
grad_u_z =@(r,z) 0;
end
