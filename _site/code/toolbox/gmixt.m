function x = gmixt(N,p,m,s)

% gmixt - Gaussian mixture
%
%   x = gmixt(N,p,s,m);
%
%   Copyright (c) 2017 Gabriel Peyre

x = randn(N,1); I = find(rand(N,1)<p); J = setdiff(1:N,I);
x(I) = x(I)*s(1)+m(1); x(J) = x(J)*s(2)+m(2); 

end