function [PC, powers] = powercentersPD(T, E, wts)
% function [PC, powers] = powercentersPD(T, E, wts)
%
% T: triangulation
% E: set of points
% wts: weights of points in E
%
% The output array PC contains the power centers of the triangles in T.
% That is, row i of PC contains the point that is equidistant from each
% vertex of the triangle specified by row i of T, with respect to the power
% of each vertex. The output array powers contains the power of each power
% center.

[m, n] = size(T);
PC = zeros(m,n-1);
powers = zeros(m,1);

for i=1:m
    tr = E(T(i,:),:);
    wt = wts(T(i,:));
    p = tr(1,:);
    Rp = repmat(p, n-1, 1);
    Pts = tr(2:n,:);
    Ac = 2*(Pts - Rp);
    
    Rw1 = repmat(wt(1), n-1, 1);
    Wts = wt(2:n);
    Sp1 = repmat(sum(p.^2), n-1, 1);
    SPts = sum(Pts.^2, 2);
    Bc = Rw1 - Wts - Sp1 + SPts;
    
    pc = Ac \ Bc;
    
    PC(i,:) = pc;
    powers(i,1) = norm(pc - p')^2 - wt(1);
end