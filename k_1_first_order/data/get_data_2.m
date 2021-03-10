function [u_vec_r,u_vec_th,u_vec_z,s,f_vec_r,f_vec_th,f_vec_z] = get_data_2(n)
%GET_DATA_2 Get u vector and s equation exact solution
%   data2
%   u = [ r^2 - r
%         0
%         z^2 - z ]
%   s = -(3r + 2z - 3)
%   f = [ -3 + (n^2/r)(r - 1)
%         (n/r)(3r - 2) - n
%         -2 - (n^2/r^2)(z^2 - z) ]
% Author: Nicole Stock
% Date: Spring 2021

u_vec_r = @(r,z) r.^2 - r;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) z.^2 - z;
s = @(r,z) -(3.*r + 2.*z - 3);
f_vec_r = @(r,z) -3 + (n.^2./r).*(r - 1);
f_vec_th = @(r,z) (n./r).*(3.*r - 2) - n;
f_vec_z = @(r,z) -2 + (n.^2./r.^2).*(z.^2 - z);
end