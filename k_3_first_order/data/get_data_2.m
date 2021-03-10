function [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_2(n)
%GET_DATA_2 Get z vector and p equation exact solution
%   data2
%   z = [ -sin(pi*z)(3r^2 - 2r)
%         -n*sin(pi*z)(r^2 - r)
%         -pi*cos(pi*z)(r^3 - r^2) ]
%   p = sin(pi*z)(r^3 - r^2)
%   f = -sin(pi*z)(6r - 2) - sin(pi*z)(3r - 2) + n^2sin(pi*z)(r - 1)
%           + pi^2sin(pi*z)(r^3 - r^2)
% Author: Nicole Stock
% Date: Spring 2021

z_vec_r = @(r,z) -sin(pi.*z).*(3.*r.^2 - 2.*r);
z_vec_th = @(r,z) -n.*sin(pi.*z).*(r.^2 - r);
z_vec_z = @(r,z) -pi*cos(pi.*z).*(r.^3 - r.^2);
p_exact = @(r,z) sin(pi.*z).*(r.^3 - r.^2);
f = @(r,z) -sin(pi.*z).*(6.*r - 2) - sin(pi.*z).*(3.*r - 2) ...
    + n.^2.*sin(pi.*z).*(r - 1) + pi.^2.*sin(pi*z).*(r.^3 - r.^2);
end