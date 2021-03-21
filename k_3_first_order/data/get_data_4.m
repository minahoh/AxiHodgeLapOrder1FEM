function [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_4(n)
%GET_DATA_4 Get z vector and p equation exact solution
%   data2
%   z = [  ]
%   p = 
%   f = 
% Author: Nicole Stock
% Date: Spring 2021

z_vec_r = @(r,z) -(z.^4 - z.^3).*(4.*r.^3 - 3.*r.^2) -sin(pi.*z).*(4.*r.^3 - 3.*r.^2);
z_vec_th = @(r,z) -n.*(z.^4 - z.^3).*(r.^3 - r.^2) -n.*sin(pi.*z).*(r.^3 - r.^2);
z_vec_z = @(r,z) -(4.*z.^3 - 3.*z.^2).*(r.^4 - r.^3)-pi.*cos(pi.*z).*(r.^4 - r.^3);
p_exact = @(r,z) (z.^4 - z.^3).*(r.^4 - r.^3) + sin(pi.*z).*(r.^4 - r.^3);
f = @(r,z) -(z.^4 - z.^3).*(12.*r.^2 - 6.*r) - (z.^4 - z.^3).*(4.*r.^2 - 3.*r) ...
    + n.^2.*(z.^4 - z.^3).*(r.^2 - r) - (12.*z.^2 - 6.*z).*(r.^4 - r.^3) ...
    -sin(pi.*z).*(12.*r.^2 - 6.*r) - sin(pi.*z).*(4.*r.^2 - 3.*r) ...
    + n.^2.*sin(pi.*z).*(r.^2 - r) + pi.^2.*sin(pi*z).*(r.^4 - r.^3);
end