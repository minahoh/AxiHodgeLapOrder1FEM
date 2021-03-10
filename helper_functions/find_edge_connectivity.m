function [t_ed] = find_edge_connectivity(t,ed,mesh)
% FIND_EDGE_CONNECTIVITY - Find the edges in each triangulation in the mesh
%
% Syntax:
%     [ed, t_ed] = find_edges(t)
% Inputs:
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumNodes matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     mesh - optional parameter indicating mesh level
%
% Outputs:
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020

if exist('mesh','var')
	% if mesh is given
    % determine if we have saved a matrix for the given mesh & domain to
    % save computational time
    if mesh == 6 && isfile('edge_resources/t_ed_6.mat')
        t_ed_6 = load('edge_resources/t_ed_6.mat');
        t_ed = t_ed_6.t_ed_6;
    elseif mesh == 7  && isfile('edge_resources/t_ed_7.mat')
        t_ed_7 = load('edge_resources/t_ed_7.mat');
        t_ed = t_ed_7.t_ed_7;
    elseif mesh == 8  && isfile('edge_resources/t_ed_8.mat')
        t_ed_8 = load('edge_resources/t_ed_8.mat');
        t_ed = t_ed_8.t_ed_8;
    elseif mesh == 9  && isfile('edge_resources/t_ed_9.mat')
        t_ed_9 = load('edge_resources/t_ed_9.mat');
        t_ed = t_ed_9.t_ed_9;
    else
        t_ed = find(t,ed);
    end
else
    t_ed = find(t,ed);
end
 
end
 
% subfunction
function [t_ed] = find(t,ed)
    [m,~] = size(ed);
    [mt,nt] = size(t);
    t_ed = NaN(mt-1,nt);

    % for each edge (rep as a row in edges_)
    for i = 1:m
        % for each triangle
        for T = 1:nt
            % if each vertex of edge is a member of triangle T
            if ismember(ed(i,1),t(1:3,T)) && ismember(ed(i,2),t(1:3,T))
                if isnan(t_ed(1,T))
                    t_ed(1,T) = i;
                elseif isnan(t_ed(2,T))
                    t_ed(2,T) = i;
                else
                    t_ed(3,T) = i;
                end                
            end
        end
    end
end