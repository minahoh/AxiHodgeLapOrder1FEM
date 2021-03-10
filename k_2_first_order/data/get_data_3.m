function [u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_3(n)
%GET_DATA_3 Get u vector and s equation exact solution
%   data3
%   s = [ n(r^3 - r^2) - (2z - 1)(r^4 - r^3)
%         -(4r^3 - 3r^2)
%         (z^2 - z)(5r^3 - 4r^2)          ]
%   u = [ 0
%         (z^2 - z)(r^4 - r^3)
%         r^4 - r^3 ]
%
% Author: Nicole Stock
% Date: Spring 2021

s_vec_r = @(r,z) n.*(r.^3 - r.^2) - (2.*z - 1).*(r.^4 - r.^3);
s_vec_th = @(r,z) -4.*r.^3 + 3.*r.^2;
s_vec_z = @(r,z) (z.^2 - z).*(r.^3 - r.^2) + (z.^2 - z).*(4.*r.^3 - 3.*r.^2);
u_vec_r = @(r,z) 0;
u_vec_th = @(r,z) (z.^2 - z).*(r.^4 - r.^3);
u_vec_z = @(r,z) r.^4 - r.^3;
f_vec_r = @(r,z) n.*(z.^2 - z).*(3.*r.^2 - 2.*r) ...
    - n.*(z.^2 - z).*(5.*r.^2 - 4.*r);
f_vec_th = @(r,z) n.^2.*(z.^2 - z).*(r.^2 - r) - 2.*(r.^4 - r.^3) ...
    - (z.^2 - z).*(15.*r.^2 - 8.*r);
f_vec_z = @(r,z) n.*(2.*z - 1).*(r.^3 - r.^2) ...
    + n.^2.*(r.^2 - r) - n.*(2.*z - 1).*(r.^3 - r.^2) ...
    - (4.*r.^2 - 3.*r) - (12.*r.^2 - 6.*r);
end