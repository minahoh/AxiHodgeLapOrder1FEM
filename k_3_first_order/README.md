# k = 3: The Dirichlet Problem for the Axisymmetric Poisson Equation

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bdiv%7D%5En_%7Brz%7D%20%5Ctext%7Bgrad%7D%5E%7Bn*%7D_%7Brz%7D%20u%20%26%3D%20f%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20u%20%26%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

## Usage

### Syntax
To compare a known exact solution u and to its approximated solution:
```
[err] = weighted_HL_k_3_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,z_vec_r,z_vec_th,z_vec_z,p_exact,n)
```
To find the approximated solution to an unknown solution:

```
[basis_p1,basis_rt1,z_h,p_h] = weighted_HL_k_3_first(f,gd,sf,ns,mesh,n)
```

### Inputs 
`f` - given function  
`gd,sf,ns` - outputs of pdepoly specifying domain  
`mesh` - max mesh level  
`z_vec_r` - exact solution z vector r component  
`z_vec_th` - exact solution z vector theta component  
`z_vec_z` - exact solution z vector z component  
`p_vec_r` - exact solution p vector r component  
`p_vec_th` - exact solution p vector theta component  
`p_vec_z` - exact solution p vector z component  

### Outputs
`err_z` - array of L2 errors for mesh levels corresponding to indices  
`err_p` - array of L2 errors for mesh levels corresponding to indices  
`basis_vertices` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_vertices(:,i,T) represents the pieceiwise basis function for the ith edge in triangle T.  
`basis_rt1` - a matrix representing piece-wise basis functions for each edge and triangle in each triangle. basis_rt1(:,i,T) represents the pieceiwise basis function for the ith edge in triangle T.  
`z_h` - approximated solution for z  
`p_h` - approximated solution for p  

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
[z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_1(n);
% run the program
[err_z,err_p] = weighted_HL_k_3_e(f,gd,sf,ns,mesh,z_vec_r,z_vec_th,z_vec_z,p_exact,n)
```

Return to [main](../README.md) page
