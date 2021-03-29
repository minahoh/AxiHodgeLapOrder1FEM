# k = 0: The Neumann Problem for the Axisymmetric Poisson Equation

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Baligned%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20-%20%5Ctext%7Bdiv%7D%5E%7Bn*%7D_%7Brz%7D%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%26%20%3D%20f%20%26%26%5Ctext%7B%20in%20%7D%20%5COmega%2C%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ctext%7Bgrad%7D%5En_%7Brz%7D%20u%20%5Ccdot%20n%20%26%20%3D%200%20%26%26%5Ctext%7B%20on%20%7D%20%5CGamma_1.%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Cend%7Baligned%7D">

## Usage

### Syntax

To compare a known exact solution u and to its approximated solution:
```
[err,grad_err,max_err] = weighted_HL_k_0_p2_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z)
```
To find the approximated solution to an unknown solution:
```
[basis,u_h] = weighted_HL_k_0_p2(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n)
```

### Inputs
`f` - given function  
`grad_f_r` - gradient(f) with respect to r  
`grad_f_z` - gradient(f) with respect to z  
`gd,sf,ns` - outputs of pdepoly specifying domain  
`mesh` - max mesh level  
`n` - n-th Fourier mode  
`u` - exact solution function  
`grad_u_r` - gradient(u) with respect to r  
`grad_u_z` - gradient(u) with respect to z  

### Outputs
`err` - array of L2 errors for mesh levels corresponding to indices  
`grad_err` - array of L2 gradient errors for mesh levels corresponding to indices  
`max_err` - array of max errors for mesh levels corresponding to indicies  
`basis` - a matrix representing piece-wise basis functions for each node in each triangle. basis(i,:,k) represents the pieceiwise basis function for the ith node in triangle k.  
`u_h` - approximated solution vector for u  

## Example
```
% add path for get_data_7() function
addpath ../data/
% define the highest mesh level
mesh = 5;
% define the nth-Fourier mode
n = 1;
% define the problem domain
pdepoly([0,1,1,0], [0,0,1,1]);
% define the equations
[f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data7();
% run the program
[err,grad_err,max_err] = weighted_HL_k_0_p2_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z);
```

Return to [main](../README.md) page
