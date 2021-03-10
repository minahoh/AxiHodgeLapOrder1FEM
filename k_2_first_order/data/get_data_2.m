function [u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_2(n)
%GET_DATA_2 Get u vector and s equation exact solution
%   data2
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
% Date: Spring 2021

s_vec_r = @(r,z) n.*(r.^3 - r.^2);
s_vec_th = @(r,z) -4.*r.^3 + 3.*r.^2;
s_vec_z = @(r,z) 0;
u_vec_r = @(r,z) 0;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) r.^4 - r.^3;
f_vec_r = @(r,z) 0;
f_vec_th = @(r,z) 0;
f_vec_z = @(r,z) n.^2.*(r.^2 - r) - 16.*r.^2 + 9.*r;
end