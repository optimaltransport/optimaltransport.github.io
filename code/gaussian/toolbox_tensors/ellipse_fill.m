function ellipse_fill(a,b,z,x0,y0,C,EC)

% ellipse_fill - draw filled ellipse
%
%   ellipse_fill(a,b,z,x0,y0,C);
%
%   Copyright (c) 2016 Gabriel Peyre

if nargin<6
    C = 'k';
end
if nargin<7
    EC = 'k';
end

q = 32;
t = linspace(0,2*pi,q);
x = a*cos(t);
y = b*sin(t);

[x,y] = deal(x0 + x*cos(z)+y*sin(z),y0 + -x*sin(z)+y*cos(z));

if not(isempty(EC))
    patch(x,y,C, 'EdgeColor', EC);
else
    patch(x,y,C);
end    
end