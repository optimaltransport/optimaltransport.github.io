%%
% Interpolation of 1D Gaussians

rep = 'results/gaussian-1d/';
[~,~] = mkdir(rep);

N = 1024;
x = linspace(0,1,N)';
normalize = @(x)x/sum(x(:));
gauss = @(m,s)normalize( exp( -(x-m).^2/(2*s^2) ) );

m0 = .1;  s0 = .025;
m1 = .8;  s1 = .08;


q = 16;
clf; hold on;
for i=1:q
    t = (i-1)/(q-1);
    m = (1-t)*m0+t*m1;
    s = (1-t)*s0+t*s1;
    lw = 2;
    plot(x, gauss(m,s), 'color', [1-t 0 t], 'LineWidth', lw);    
end
set(gca, 'XTick', [], 'YTick', []);
box on; axis tight;
saveas(gcf, [rep 'interp-density.eps']);


clf; hold on;
for i=1:q
    t = (i-1)/(q-1);
    m = (1-t)*m0+t*m1;
    s = (1-t)*s0+t*s1;
    lw = 2;
    plot(x, gauss(m,s), 'color', [1-t 0 t], 'LineWidth', lw);    
end
set(gca, 'XTick', [], 'YTick', []);
box on; axis tight;
saveas(gcf, [rep 'interp-density.eps'], 'epsc');