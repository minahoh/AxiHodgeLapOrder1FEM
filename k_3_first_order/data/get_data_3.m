function [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_3(n)
%GET_DATA_2 Get z vector and p equation exact solution
%   data2
%   z = [ -sin(pi*z)(4r^3 - 3r^2)
%         -n*sin(pi*z)(r^3 - r^2)
%         -pi*cos(pi*z)(r^4 - r^3) ]
%   p = sin(pi*z)(r^4 - r^3)
%   f = -sin(pi*z)(12r^2 - 6r) - sin(pi*z)(4r^2 - 3r) 
%           + n^2sin(pi*z)(r^2 - r) + pi^2sin(pi*z)(r^4 - r^3)
% Author: Nicole Stock
% Date: Spring 2021

z_vec_r = @(r,z) -sin(pi.*z).*(4.*r.^3 - 3.*r.^2);
z_vec_th = @(r,z) -n.*sin(pi.*z).*(r.^3 - r.^2);
z_vec_z = @(r,z) -pi.*cos(pi.*z).*(r.^4 - r.^3);
p_exact = @(r,z) sin(pi.*z).*(r.^4 - r.^3);
f = @(r,z) -sin(pi.*z).*(12.*r.^2 - 6.*r) - sin(pi.*z).*(4.*r.^2 - 3.*r) ...
    + n.^2.*sin(pi.*z).*(r.^2 - r) + pi.^2.*sin(pi*z).*(r.^4 - r.^3);
end