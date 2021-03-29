function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data8()
%GET_U_DATA8 Get u and f function handles
%   data3
%   u = r^(4/3)
% Author: Nicole Stock
% Date: Fall 2020

u = @(r,z) r.^(4./3);
grad_u_r =@(r,z) (4./3).*r.^(1./3);
grad_u_z =@(r,z) 0;

f = u;
grad_f_r = grad_u_r;
grad_f_z = grad_u_z;
end
