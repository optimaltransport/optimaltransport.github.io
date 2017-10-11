function [vxx,vy] = powerdiagram(varargin)
% Power diagram
% modified voronoi.m to use power centers instead of circumcenters

[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,5,nargs));

x = args{1};
y = args{2};
if ~isequal(size(x),size(y))
        error(message('MATLAB:voronoi:InputSizeMismatch'));
end
if ndims(x) > 2 || ndims(y) > 2
        error(message('MATLAB:voronoi:HigherDimArray'));
end   
x = x(:);
y = y(:);
tri = args{3};
ls = args{4};
wts = args{5};

if isempty(tri)
    return;
end

% Compute power centers of triangles
% (this is already done in cells2D)
sE = [x, y];
[c, powers] = powercentersPD(tri, sE, wts);
% tr = triangulation(tri,x,y);
% c = tr.circumcenter();


% Create matrix T where i and j are endpoints of edge of triangle T(i,j)
n = numel(x);
t = repmat((1:size(tri,1))',1,3);
T = sparse(tri,tri(:,[3 1 2]),t,n,n); 

% i and j are endpoints of internal edge in triangle E(i,j)
E = (T & T').*T; 
% i and j are endpoints of external edge in triangle F(i,j)
F = xor(T, T').*T;

% v and vv are triangles that share an edge
[~,~,v] = find(triu(E));
[~,~,vv] = find(triu(E'));

% Internal edges
vx = [c(v,1) c(vv,1)]';
vy = [c(v,2) c(vv,2)]';

%%% Compute lines-to-infinity
% i and j are endpoints of the edges of triangles in z
[i,j,z] = find(F);
% Counter-clockwise components of lines between endpoints
dx = x(j) - x(i);
dy = y(j) - y(i);

% Calculate scaling factor for length of line-to-infinity
% Distance across range of data
rx = max(x)-min(x); 
ry = max(y)-min(y);
% Distance from vertex to center of data
cx = (max(x)+min(x))/2 - c(z,1); 
cy = (max(y)+min(y))/2 - c(z,2);
% Sum of these two distances
nm = sqrt(rx.*rx + ry.*ry) + sqrt(cx.*cx + cy.*cy);
% Compute scaling factor
scale = nm./sqrt((dx.*dx+dy.*dy));
    
% Lines from voronoi vertex to "infinite" endpoint
% We know it's in correct direction because compononents are CCW
ex = [c(z,1) c(z,1)-dy.*scale]';
ey = [c(z,2) c(z,2)+dx.*scale]';
% Combine with internal edges
vx = [vx ex];
vy = [vy ey];

if nargout<2
    % Plot diagram
    if isempty(cax)
        % If no current axes, create one
        cax = gca;
    end
    if isempty(ls)
        % Default linespec
        ls = '-';
    end
    [l,c,mp,msg] = colstyle(ls); error(msg) % Extract from linespec
    if isempty(mp)
        % Default markers at points        
        mp = '.';
    end
     if isempty(l)
        % Default linestyle
        l = get(ancestor(cax,'figure'),'DefaultAxesLineStyleOrder'); 
    end
    if isempty(c), 
        % Default color        
        co = get(ancestor(cax,'figure'),'DefaultAxesColorOrder');
        c = co(1,:);
    end
    % Plot points
    h1 = plot(x,y,'marker',mp,'color',c,'linestyle','none','parent',cax);
    % Plot voronoi lines
    h2 = line(vx,vy,'color',c,'linestyle',l,'parent',cax,...
        'yliminclude','off','xliminclude','off');
    if nargout==1, vxx = [h1; h2]; end % Return handles
else
    vxx = vx; % Don't plot, just return vertices
end