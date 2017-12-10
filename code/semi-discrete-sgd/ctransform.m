function [v,g] = ctransform(x, u, y, p, epsilon)

% ctransform - c-transform and gradient of semi-discrete energy
%
%   [v,g] = ctransform(x,u, y, epsilon);
%
%   u is sampled at points x. 
%   v is u^c sampled at points y. 
%
%   x must be of size (dim,n)
%   y must be of size (dim,m)
%
% If epsilon==0
%   u^c(y) = min_x  |x-y|^p-u(x)
% else
%   u^c(y) = - epsilon * log( sum_x  exp( (-|x-y|^p + u(x))/epsilon )
%
%   The semi discrete energy between continuous measure a and discrete
%   measure b is
%       E(u) = int u^c(y) da(y) - <u,b>
%   g(:,j) is the gradient of u -> u^c(y_j).
%
%   Copyright (c) 2017 Gabriel Peyre

if nargin<5
    epsilon = 0;
end

n = size(x,2);
m = size(y,2);

c = distmat(x,y).^p; % cost
S = c -  repmat(u(:), [1 m]);
[v,I] = min(S);
if epsilon>0
    % use soft-min, stabilized
    if 1
        S = S - repmat(v, [size(S,1) 1]);
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) ) + v;
    else
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) );    
    end
end
v = v(:);

if nargout>1
    % compute also the gradient
    if epsilon==0
        % simply indicator function of the Laguerre cells, indicated by the index I
        g = sparse( I, 1:m, ones(m,1),n, m );
    else
        % smoothed-indicator
        g = exp( -S/epsilon );
        g = g ./ repmat( sum(g), [n 1] );
    end
    
end

end
