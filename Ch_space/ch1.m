function [err] = ch1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
%CH1 Fourier Raviart Thomas and P1 Space Finite 
%  Element Method. 
%   This program is set up to be given an exact solution.
%   First Order Fourier Raviart Thomas and P1 (Ch_1) Space
%
% Syntax:
%     [err] = ch1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
%     f_vec_r - given vector r component
%     f_vec_z - given vector z component
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     u_vec_r - exact solution vector r component
%     u_vec_z - exact solution vector z component
%
% Outputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    addpath ../helper_functions data
%    n = 1;
%    [u_vec_r,u_vec_th,u_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err] = ch1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n);
% Dependencies:
%    ../new_ele[mesh].mat
%    basis_functions_weighted_HL_k_0_p1.m
%    basis_functions_rt1.m
%    create_b_ch_1.m
%    display_errors.m
%    errors_exact_ch_1.m
%    stiffness_matrix_ch_1.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
[~,triangles]=size(t);
ele=t(1:3,1:triangles);
ele=ele';
node=p';
tr=triangulation(ele,node);
ed=edges(tr);
load(['../new_ele',num2str(1),'.mat']);

% Init error vector
err = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   

    [basis_vertices,basis_rt1,x] = solve(p,t,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n);
    err(1) = errors_exact_ch_1(p,t,ed,new_ele,basis_vertices,basis_rt1,x,u_vec_r,u_vec_th,u_vec_z,n);
    
    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        [~,triangles]=size(t);
        ele=t(1:3,1:triangles);
        ele=ele';
        node=p';
        tr=triangulation(ele,node);
        ed = edges(tr);
        load(['../new_ele',num2str(i),'.mat']);
        
        [basis_vertices,basis_rt1,x] = solve(p,t,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n);
        err(i) = errors_exact_ch_1(p,t,ed,new_ele,basis_vertices,basis_rt1,x,u_vec_r,u_vec_th,u_vec_z,n);
    end
    display_errors(err);

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis_vertices,basis_rt1,x] = solve(p,t,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n)
    basis_vertices = basis_functions_weighted_HL_k_0_p1(p,t);
    basis_rt1 = basis_functions_rt1(p,t,ele,ed,new_ele);
    S = stiffness_matrix_ch_1(p,t,ed,new_ele,basis_vertices,basis_rt1,n);
    b = create_b_ch_1(p,t,ed,new_ele,basis_vertices,basis_rt1,f_vec_r,f_vec_th,f_vec_z,n);
    x = S\b;

end
