function v = ctransform(c, x,u,y, epsilon)

% ctransform - c-transform 
%
%   v = ctransform(c, x,u,y, epsilon);
%
% If epsilon==0
%   u^c(y) = min_x  c(x,y)-u(x)
% else
%   u^c(y) = - epsilon * log( sum_x  exp( (-c(x,y)+u(x))/epsilon )
%
%   Copyright (c) 2017 Gabriel Peyre

if nargin<5
    epsilon = 0;
end


[Y,X] = meshgrid(y,x);
S = c(X,Y) -  repmat(u(:), [1 length(y)]);
v = min(S);
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

end
