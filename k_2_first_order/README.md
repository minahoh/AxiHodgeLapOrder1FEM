# k = 2: The Axisymmetric Vector Laplacian curl curl + grad div

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bcurl%7D%5En_%7Brz%7D%20%5Ctext%7Bcurl%7D%5E%7Bn*%7D_%7Brz%7D%20u%20-%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%20f%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7Brz%7D%20%5Ccdot%20t%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u_%7B%5Ctheta%7D%20%26%3D%200%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20u%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

## Usage


### Syntax
To compare a known exact solution u and to its approximated solution:
```
[err] = weighted_HL_k_2_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,n)
```
To find the approximated solution to an unknown solution:
```
[basis_p1,basis_rt1,basis_p2,basis_nd1,u_h,s_h] = weighted_HL_k_2_first(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n)
```

### Input

`f_vec_r` - given function r component  
`f_vec_th` - given function theta component  
`f_vec_z` - given function z component  
`gd,sf,ns` - outputs of pdepoly specifying domain  
`mesh` - max mesh level  
`u_vec_r` - exact solution z vector r component  
`u_vec_th` - exact solution z vector theta component  
`u_vec_z` - exact solution z vector z component  
`s_vec_r` - exact solution s vector r component  
`s_vec_r` - exact solution s vector theta component  
`s_vec_z` - exact solution s vector z component  

### Outputs
`err_u` - array of L2 errors for mesh levels corresponding to indices  
`err_s` - array of L2 errors for mesh levels corresponding to indices   
`basis_p1` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_vertices(:,i,T) represents the pieceiwise basis function for the ith edge in triangle T.  
`basis_rt1` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_rt1(:,i,T) represents the pieceiwise basis function for the ith edge in triangle T.  
`basis_p2` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_p2(:,i,T) represents the pieceiwise basis function for the ith node or midpoint in triangle T.  
`basis_nd1` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_nd1(:,i,T) represents the ith pieceiwise basis function in triangle T.  
`u_h` - approximated solution for u  
`s_h` - approximated solution for s  

## Example
```
% add path for get_data_1() function
addpath ../data/
% define the highest mesh level
mesh = 5;
% define the nth-Fourier mode
n = 1;
% define the problem domain
pdepoly([0,1,1,0], [0,0,1,1]);
% define the equations
[u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
% run the program
[err_u,err_s] = weighted_HL_k_2_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,n);
```

Return to [main](../README.md) page
