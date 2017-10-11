function L = sinkhorn_loss(x,y,mu,nu,epsilon,p, sinkh_mode, options)

% sinkhorn_loss - sinkhorn loss between two discrete measures
%
%   L = sinkhorn_loss(x,y,mu,nu,epsilon,p, sinkh_mode, options)
%
%   sinkh_mode=='uncorrected' (default)
%   sinkh_mode=='mmd': apply a correction to be consisten with Energy
%   Distance when epsilon->+inf
% 
%   Copyright (c) 2017 Gabriel Peyre

n = size(x,2);

options.null = 0;
c = distmat(x,y).^p;
options.verb = 0;
options.tol = 1e-10;
[u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
L = sum(gamma(:).*c(:));

if strcmp(sinkh_mode, 'mmd');
    %% Self costs %%
    c = distmat(x,x).^p;
    [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,mu,c,epsilon,options);
    Lx = sum(gamma(:).*c(:));
    %
    c = distmat(y,y).^p;
    [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(nu,nu,c,epsilon,options);
    Ly = sum(gamma(:).*c(:));
    L = L - Lx/2 - Ly/2;
end

end