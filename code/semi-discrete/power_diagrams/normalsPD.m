function ind = normalsPD(LE, C)
% function ind = normalsPD(E, C)
%
% LE: set of lifted points
% C: convex hull of LE
%
% The output ind contains indices of facets of the lower hull.

[m, n] = size(C);
center = mean(LE,1);

for i=1:m
    v = null(bsxfun(@minus, LE(C(i,1),:), LE(C(i,2:end),:)))';
    [mm, nn] = size(v);
    if mm > 1
        V(i,:) = NaN;
        error('nullspace error')
        % possibility of degenerate null vectors
    else
        V(i,:) = v;
    end
    mid(i,:) = mean(LE(C(i,:),:),1);
end

dot = sum(bsxfun(@minus, center, mid).*V, 2);
outer = dot < 0;
V(outer,:) = -1*V(outer,:);

ind = V(:,n) > 0;