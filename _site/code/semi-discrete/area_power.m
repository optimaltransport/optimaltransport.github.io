function A = area_power(xy,V,C,bb)

% area_power - area of a power diagram

N = size(xy,1);

% plot the sites' points
A = [];
for i = [1:N]
    if not(isempty(C{i}))
        xpoly = V(C{i},1);
        ypoly = V(C{i},2);
        xpoly = [xpoly; xpoly(1)];
        ypoly = [ypoly; ypoly(1)];
        A(i) = polyarea(xpoly,ypoly);
    else
        A(i) = 0;
    end
end

end