function [u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n)
%GET_DATA1 Get u vector and s equation exact solution
%   data1
%   s = [ nr(r-1)
%         -3r^2 + 2r
%         0          ]
%   u = [ 0
%         0
%         r^2(r-1) ]
%   f = [ 0
%         0
%         rn^2 - n^2 - 9r + 4 ]
% Author: Nicole Stock
% Date: Fall 2020

s_vec_r = @(r,z) n.*r.*(r-1);
s_vec_th = @(r,z) -3.*r.^2 + 2.*r;
s_vec_z = @(r,z) 0;
u_vec_r = @(r,z) 0;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) (r.^2).*(r-1);
f_vec_r = @(r,z) 0;
f_vec_th = @(r,z) 0;
f_vec_z = @(r,z) (n.^2).*r - (n.^2) - 9.*r + 4;
end