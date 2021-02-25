function [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_1(n)
%GET_DATA1 Get z vector and p equation exact solution
%   data1
%   z = [ -sin(pi*z)(2r - 1)
%         -sin(pi*z)(r - 1)
%         -pi*cos(pi*z)(r^2 - r) ]
%   p = sin(pi*z)(r^2 - r)
%   f = -2sin(pi*z) - (1/r)sin(pi*z)(2r - 1 - n^2(r - 1)) -
%       pi^2cos(pi*z)(r^2 - r)
% Author: Nicole Stock
% Date: Fall 2020

z_vec_r = @(r,z) -sin(pi.*z).*(2.*r - 1);
z_vec_th = @(r,z) -sin(pi.*z).*(r - 1);
z_vec_z = @(r,z) -pi*cos(pi.*z).*(r.^2 - r);
p_exact = @(r,z) sin(pi.*z).*(r.^2 - r);
f = @(r,z) -2.*sin(pi.*z) ...
    - (1./r).*sin(pi.*z).*(2.*r - 1 - n.^2.*(r - 1)) ...
    + (pi.^2).*sin(pi.*z).*(r.^2 - r);
end