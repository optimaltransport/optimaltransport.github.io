function d = rkhs_norm_invariant(k,x,p,y,q)

% compute a dual of RKHS norm on measures
%
%   d = rkhs_norm_invariant(k,x,p,y,q);
%
%   d is the *squared* RKHS dual norm
%   d = sum_ii' pi*pi'*k(|xi-xi'|) + sum_jj' qj*qj'*k(|xj-xj'|) - 2 sum_ij pi*qj*k(|xi-yj|)
%
%   Copyright (c) 2017 Gabriel Peyre

p = p(:); q = q(:);

D = distmat(x,x);
A = (p*p').*k(D);
%
D = distmat(y,y);
B = (q*q').*k(D);
%
D = distmat(x,y);
C = (p*q').*k(D);

d = sum(A(:))+sum(B(:))-2*sum(C(:));

end
