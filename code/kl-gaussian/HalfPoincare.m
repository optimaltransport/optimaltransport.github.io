%%
% plot geodesic on the Poincar? half plane

rep = '../results/kl-gaussian/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);


%%
% Distance on Poincarr? space.

D = @(u,v)acosh(1 + abs(u-v).^2 ./ (2*imag(u).*imag(v))  );

u = 0 + .5i;

m = linspace(-1,1,200);
s = linspace(1e-3,2,200);
[S,M] = meshgrid(s,m);
V = M + 1i*S;

clf; hold on;
imagesc(m,s,D(u,V)');
contour(m,s,D(u,V)', 12, 'k');
colormap parula(256);
axis xy;


%%
% Generate circle

r = 1000;
tmin = .3; tmax = .95;
tmin = .15; tmax = .98;
t = linspace(tmin*pi,tmax*pi,r);
v = exp(1i*t);

% re-sample at equi-spacing
s = D(v(1),v);
s = s/s(end);
q = 15; % #intermediate points
ti = interp1(s, linspace(0,1,r), linspace(0,1,q));
ti = tmin + (tmax-tmin)*ti;
vi = exp(1i*pi*ti);


clf; hold on;
%plot( exp(1i*linspace(0,pi,r)) , 'k:');
plot(v, 'color', 'k', 'LineWidth', 2);
for i=1:q
    t = (i-1)/(q-1);
    col = [t 0 1-t];
    plot(vi(i), '.', 'MarkerSize', 45, 'MarkerEdgeColor', col);
end
axis equal; axis([-1 1 0 1]); 
axis off;
saveas(gcf, [rep 'geod-kl.eps'], 'epsc');


t = linspace(0,1,r);
t1 = linspace(-.5,1.5,r);
clf; hold on;
plot( v(1)*t+v(end)*(1-t) , 'color', 'k', 'LineWidth', 2);
%plot( v(1)*t1+v(end)*(1-t1) , 'k:');
for i=1:q
    t = (i-1)/(q-1);
    col = [t 0 1-t];
    plot(v(1)*(1-t)+v(end)*t, '.', 'MarkerSize', 45, 'MarkerEdgeColor', col);
end
axis equal; axis([-1 1 0 1]); 
axis off;
saveas(gcf, [rep 'geod-ot.eps'], 'epsc');


%%
% interpolations of densities. 

x = linspace(-2,2.5, 1024);
gauss = @(m,s)1/s * exp(-(x-m).^2/(2*s^2));

clf; hold on;
for i=1:q
    t = (i-1)/(q-1);
    vt = v(1)*(1-t)+v(end)*t;
    m = sqrt(2)*real(vt);
    s = imag(vt);
    col = [t 0 1-t];
    plot(x, gauss(m,s), 'Color', col, 'LineWidth', 2);
end
axis tight; axis off;
SetAR(1/2);
saveas(gcf, [rep 'densities-ot.eps'], 'epsc');

clf; hold on;
for i=1:q
    t = (i-1)/(q-1);
    vt = v(1)*(1-t)+v(end)*t;
    m = sqrt(2)*real(vi(i));
    s = imag(vi(i));
    col = [t 0 1-t];
    plot(x, gauss(m,s), 'Color', col, 'LineWidth', 2);
end
axis tight; axis off;
SetAR(1/2);
saveas(gcf, [rep 'densities-kl.eps'], 'epsc');

