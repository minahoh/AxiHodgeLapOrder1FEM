function [u_vec_r,u_vec_th,u_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n)
%GET_DATA_1 Get u vector and f vector
%   u = [ z - (1/n)((1/3)r^3 - (1/2)r^2) ]
%       [  -nz + (1/3)r^3 - (1/2)r^2     ]
%       [ r                              ]
%   F = [ n(r-1) + u_r   ]
%       [ -2r + 1 + u_th ]
%       [ u_z            ]

u_vec_r = @(r,z) z - (1./n).*((1./3).*r.^3 - (1./2).*r.^2);
u_vec_th = @(r,z) (-1).*n.*z + (1./3).*r.^3 - (1./2).*r.^2;
u_vec_z = @(r,z) r;

f_vec_r = @(r,z) n.*(r-1) + z - (1./n).*((1./3).*r.^3 - (1./2).*r.^2);
f_vec_th = @(r,z) -2.*r + 1 - n.*z + (1./3).*r.^3 - (1./2).*r.^2;
f_vec_z = @(r,z) r;
end