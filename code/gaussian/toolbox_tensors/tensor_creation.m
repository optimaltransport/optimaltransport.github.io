function g = tensor_creation(r, theta)

% tensor_creation - generate tensor list from polar coordinates
%
%   g = tensor_creation(r, theta)
%
%   r in [0,1] x [0,2*pi[
%
%   a = r(i)*cos(theta(i))+1)/2;
%   b = r(i)*sin(theta(i))/2;
%   g(i,:,:) = [a b; b 1-a];
%
%   Copyright (c) Gabriel Peyre

% r = r(:); theta = theta(:);
if min(r(:))<0 || max(r(:))>1 || min(theta(:))<0 || max(theta(:))>2*pi
%    warning('Tensor parameterization error');
end

g = zeros(2,2,size(theta,1),size(theta,2));

a = (r.*cos(theta)+1)/2;
b = r.*sin(theta)/2;
resh = @(x)reshape(x, [1 1 size(x,1) size(x,2)]);
a = resh(a); b = resh(b);

g(1,1,:,:) = a; 
g(2,2,:,:) = 1-a;
g(1,2,:,:) = b; 
g(2,1,:,:) = b; 

end