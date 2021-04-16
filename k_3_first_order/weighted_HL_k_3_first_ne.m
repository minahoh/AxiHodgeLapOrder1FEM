function [err_p] = weighted_HL_k_3_first_ne(f,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_3_FIRST_NE Hodge Laplacian k = 3 First Order Finite Element Method.
%   This program is set up to be given an exact solution.
%   Hodge Laplacian k = 3 case, first order
%   {psi_i}i=1->2Ne+2Nt+N is the basis for Ch1
%   {chi_j}j=1->2Nt is the basis for Dh1
%   (Ch is the weighted fourier Raviart Thomas space)
%   (Dh is the piecewise constant space)
%   Solve for (z,p) in (Ch x Dh) s.t.
%       (z , w)_r - (p , div_rz^n(w))_r = 0
%       (div_rz^n(z) , s)_r = (f , s)_r
%           for all w in Ch, s in Dh
%
% Syntax:
%     [err] = weighted_HL_k_3_first_ne(f,gd,sf,ns,mesh,n)
%     f - given function
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     err_z - array of L2 errors for mesh levels corresponding to indices
%     err_p - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    addpath ../helper_functions data
%    n = 1;
%    [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err_p] = weighted_HL_k_3_first_ne(f,gd,sf,ns,mesh,n);
% Dependencies:
%    ../new_ele[mesh].mat.m
%    basis_functions_HL_k_3_first.m
%    create_B_HL_k_3_first.m
%    create_F_HL_k_3_first.m
%    display_errors.m
%    errors_exact_HL_k_3_first.m
%    mass_matrix_HL_k_3_first.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../');
addpath('../LShape_new_ele_files/');

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf,'MesherVersion','R2013a');
[~,triangles]=size(t);
ele=t(1:3,1:triangles);
ele=ele';
node=p';
tr=triangulation(ele,node);
ed=edges(tr);
%load(['LShapeDomain/new_ele',num2str(1),'.mat']);
%t_ed = new_ele;
%load(['t_ed_',num2str(1),'.mat']);
%  t_ed = find_edge_connectivity(t,ed);
%  t_ed = t_ed';
aa = load(['new_ele',num2str(1),'.mat']);
t_ed = cell2mat(struct2cell(aa));

pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

% Init error vector
err_p = zeros(1,mesh);
grad_err = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end

    [basis_p1,basis_rt1,z_h,p_h] = solve(p,t,ele,ed,t_ed,f,n);
    
    for i = 2:mesh    
        disp(i)
        p_h_km1 = p_h;
        
        % Refine mesh to next level
        [p,e,t,p_h_km1]=refinemesh(g,p,e,t,p_h_km1,it,'regular');
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [~,triangles]=size(t);
        ele=t(1:3,1:triangles);
        ele=ele';
        node=p';
        tr=triangulation(ele,node);
        ed = edges(tr);
        %load(['../new_ele',num2str(i),'.mat']);
        %load(['LShapeDomain/new_ele',num2str(i),'.mat']);
        %t_ed = new_ele;
        %disp(new_ele)
        %load(['t_ed_',num2str(i),'.mat']);
%          t_ed = find_edge_connectivity(t,ed);
%          t_ed = t_ed';
        aa = load(['new_ele',num2str(i),'.mat']);
        t_ed = cell2mat(struct2cell(aa));

        [basis_p1,basis_rt1,z_h,p_h] = solve(p,t,ele,ed,t_ed,f,n);
        [err_p(i),grad_err(i)] = errors_no_exact_HL_k_3_first(p,t,basis_p1,p_h,p_h_km1);
    end
    fprintf('p\n');
    display_errors(err_p,grad_err);

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis_p1,basis_rt1,z_h,p_h] = solve(p,t,ele,ed,new_ele,f,n)
    [basis_p1,basis_rt1] = basis_functions_HL_k_3_first(p,t,ele,ed,new_ele);
    M = mass_matrix_HL_k_3_first(p,t,ed,new_ele,basis_p1,basis_rt1,n);
    B = create_B_HL_k_3_first(p,t,ed,new_ele,basis_p1,basis_rt1);
    F = create_F_HL_k_3_first(p,t,basis_p1,f);
    
    Bt = B.';
    MB = (M\B);
    p_h = (Bt*MB)\F;
    z_h = MB*p_h;
end