function [u_vec_r,u_vec_th,u_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_2(n)
%GET_DATA_2 Get u vector and f vector
%   u = [ (1/2)rz - (1/2)rz^2
%         0
%         (1/3)rz^3 - (1/2)rz^2 ]
%   F = [ (1/2)rz - (1/2)rz^2 - z^2 + z
%         -(n/r)(rz^2 - z^2 - rz + z)
%         (1/3)rz^3 - (1/2)rz^2 -2rz + 2z + r - 1 ]
% Author: Nicole Stock
% Date: Fall 2020

u_vec_r = @(r,z) (1./2).*r.*z - (1./2).*r.*z.^2;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) (1./3).*r.*z.^3 - (1./2).*r.*z.^2;

f_vec_r = @(r,z) (1./2).*r.*z - (1./2).*r.*z.^2 - z.^2 + z;
f_vec_th = @(r,z) -(n./r).*(r.*z.^2 - z.^2 - r.*z + z);
f_vec_z = @(r,z) (1./3).*r.*z.^3 - (1./2).*r.*z.^2 - 2.*r.*z + 2.*z + r - 1;
end