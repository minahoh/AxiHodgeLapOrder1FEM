function [t_ed] = get_edge_connectivity(gd,sf,ns,mesh)
% GET_EDGE_CONNECTIVITY  build edge connectivitity matrix for given mesh level
%
% Syntax:
%     [t_ed] = get_edges(gd,sf,ns,mesh)
%
% Inputs:
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%
% Outputs:
%    t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%
% Usage Example:
%    mesh = 8;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [t_ed] = get_edge_connectivity(gd,sf,ns,mesh)
%
% Dependencies:
%    find_edges.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('../edge_resources')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
tr = triangulation(t(1:3,:)',p');

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   
    
    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
    end
    
    % get edge information for highest mesh level
    ed = edges(tr);
    t_ed = find_edge_connectivity(t,ed);
end

% mesh level must be greater than 1

% end main
end

