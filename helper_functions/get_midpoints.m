function [p2,t2] = get_midpoints(gd,sf,ns,mesh)
%GET_MIDPOINTS
%
% Syntax:
%     [p2,t2] = get_midpoints(gd,sf,ns,mesh)
%
% Inputs:
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%
% Outputs:
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%
% Usage Exampled:
%    mesh = 8;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [p2,t2] = get_midpoints(gd,sf,ns,mesh)
%
% Dependencies:
%    find_midpoints.m
%
% Author: Nicole Stock
% Date: Fall 2020

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
%pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

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
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');
    end
    
    % Find the midpoints for P2 nodal points for the highest mesh level
    [p2,t2] = find_midpoints(p,t);
end
% mesh level must be greater than 1

% end main
end
