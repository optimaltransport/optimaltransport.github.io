function EA = edgeAttPD(T, edges)
% function EA = edgeAttPD(T, edges)
%
% T: pieces of a triangulation
% edges: array of edges in the triangulation
% 
% Each entry in the output cell EA corresponds to an edge from the input
% edges and contains all pieces of the triangulation attached to that edge.

[me, ne] = size(edges);
EA = cell(me,1);

for i=1:me
    EA{i} = find(sum(ismember(T, edges(i,:)), 2) == ne)';
end