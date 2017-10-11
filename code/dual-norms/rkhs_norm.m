function d = rkhs_norm(k,x,p,y,q)

% compute a dual of RKHS norml on measures
%
%   d = rkhs_norm(x,p,y,q);
%
%   d is the *squared* RKHS dual norm
%   d = sum_ii' pi*pi'*k(xi,xi') + sum_jj' qj*qj'*k(xj,xj') - 2 sum_ij pi*qj*k(xi,yj)
%
%   Copyright (c) 2017 Gabriel Peyre

nx = length(p);
ny = length(q);


[J,I] = meshgrid(1:nx,1:nx);
A = p(I).*p(J).*k(x(I),x(J));
%
[J,I] = meshgrid(1:ny,1:ny);
B = q(I).*q(J).*k(y(I),y(J));
%
[J,I] = meshgrid(1:ny,1:nx);
r = @(z)reshape(z, [nx ny]);
C = r(p(I)) .* r(q(J)) .*k( r(x(I)), r(y(J)) );

d = sum(A(:))+sum(B(:))-2*sum(C(:));


end