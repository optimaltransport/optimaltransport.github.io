function [B,err] = gbary(C,lambda,niter)

% gbary - compute Wassertein barycenter of Gaussian matrices.
%
% [B,err] = gbary(C,lambda,niter);
%
%   Copyright (c) 2017 Gabriel Peyre

d = size(C{1},1);
K = length(C);
if nargin<3
    niter = 100;
end

B = eye(d);
err = [];
for it=1:niter
    U = zeros(d);
    for k=1:K
        U = U + lambda(k)*sqrtm( sqrtm(B)*C{k}*sqrtm(B) );
    end
    err(it) = norm( B-U,'fro' );
    B = U;
end

end