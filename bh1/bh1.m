function [err] = bh1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
%BH1 Fourier Raviart Thomas and P1 Space Finite 
%  Element Method. 
%   This program is set up to be given an exact solution.
%   First Order Fourier Nedelec and P2 (Bh_1) Space
%
% Syntax:
%     [err] = bh1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
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
%    [err] = bh1(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n);
% Dependencies:
%    ../new_ele[mesh].mat
%    basis_functions_weighted_HL_k_0_p2.m
%    basis_functions_bh1.m
%    create_b_bh_1.m
%    display_errors.m
%    errors_exact_bh_1.m
%    stiffness_matrix_bh_1.m
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

% Find the midpoints for P2 nodal points
[p2,t2] = find_midpoints(p,t);

% Init error vector
err = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   

    [basis_p2,basis_nd1,x] = solve(p,t,p2,t2,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n);
    err(1) = errors_exact_bh_1(p,t,t2,ed,new_ele,basis_p2,basis_nd1,x,u_vec_r,u_vec_th,u_vec_z,n);
    
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
        
        % Find the midpoints for P2 nodal points
        [p2,t2] = find_midpoints(p,t);
        
        [basis_p2,basis_nd1,x] = solve(p,t,p2,t2,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n);
        err(i) = errors_exact_bh_1(p,t,t2,ed,new_ele,basis_p2,basis_nd1,x,u_vec_r,u_vec_th,u_vec_z,n);
    end
    display_errors(err);

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis_p2,basis_nd1,x] = solve(p,t,p2,t2,ele,ed,new_ele,f_vec_r,f_vec_th,f_vec_z,n)
    basis_p2 = basis_functions_weighted_HL_k_0_p2(p,t,p2,t2);
    basis_nd1 = basis_functions_nd1(p,t,ele,ed,new_ele);
    S = stiffness_matrix_bh_1(p,t,p2,t2,ed,new_ele,basis_p2,basis_nd1,n);
    b = create_b_bh_1(p,t,p2,t2,ed,new_ele,basis_p2,basis_nd1,f_vec_r,f_vec_th,f_vec_z,n);
    x = S\b;
end
