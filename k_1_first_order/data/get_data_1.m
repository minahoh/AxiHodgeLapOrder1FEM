function [u_vec_r,u_vec_th,u_vec_z,s,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n)
%GET_DATA1 Get u vector and s equation exact solution
%   data1
%   u = [ r^3(r-1)  = r^4 - r^3
%         0
%         0       ]
%   s = -5r^3 + 4r^2
%   f = [ -15r^2 + 8r + n^2r^2 - n^2r
%         2nr^2 - 2nr
%         0                           ]
% Author: Nicole Stock
% Date: Fall 2020

u_vec_r = @(r,z) r.^4 - r.^3;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) 0;
s = @(r,z) -5.*r.^3 + 4.*r.^2;
f_vec_r = @(r,z) -15.*r.^2 + 8.*r + (n.^2).*r.^2 - (n.^2).*r;
f_vec_th = @(r,z) 2.*n.*r.^2 - 2.*n.*r;
f_vec_z = @(r,z) 0;
end