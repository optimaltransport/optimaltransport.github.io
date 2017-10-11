function disp_distrib(xi,m)

xi = xi/max(abs(xi(:)));


if nargin==1
m = 20;
end

clf; hold on;
imagesc(xi);
if m>0
    contour(xi, linspace(-1,1,m), 'k');
end
axis image; axis off;
caxis([-1 1]);
t = linspace(0,1,128)';
colormap( [ [t*0+1,t,t] ; [1-t,1-t,t*0+1] ] );

end