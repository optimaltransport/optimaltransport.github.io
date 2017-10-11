%%
% Test for c-transform in 1D

addpath('../toolbox/');

rep = 'results/c-transform/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);


kx = 10; % #samples
ky = 12;
n = 512/2; % for dense display.

randn('state', 66);
rand('state', 66);

x = rescale( cumsum( rand(kx,1)+.5 ), .05,.95 );
y = rescale( cumsum( rand(ky,1)+.5 ), .05,.95 );
g = linspace(0,1,n);
gx = sort(union(g,x));
gy = sort(union(g,y));

u = rand(kx,1);
u(end) = 1;
u(end/2) = 1;
u(1) = .9;

%% 
% Just checking

r = 1;
c = @(x,y)5*abs(x-y).^r;
epsilon = 0;
%
v = ctransform(c, x,u,y, epsilon);
vg = ctransform(c, x,u,gy, epsilon);
%
uu = ctransform(c, y,v,x, epsilon);
uug = ctransform(c, y,v,gx, epsilon);
%
vv = ctransform(c, x,uu,y, epsilon);
vvg = ctransform(c, x,uu,gy, epsilon); % should match vg



rlist = [.5 1 1.5 2];
rlist = [1];

m = .2;

eps_list = linspace(0,.5,7);

for i=1:length(rlist)
    r = rlist(i);
    c = @(x,y)5*abs(x-y).^r;

    clf; hold on;
    for ieps=1:length(eps_list)
        epsilon = eps_list(ieps);
        mc = (ieps-1)/(length(eps_list)-1);
        col = [0 mc 1-mc];
        %
        v = ctransform(c, x,u,y, epsilon);
        vg = ctransform(c, x,u,gy, epsilon);
        %
        plot(gy,vg, '-', 'LineWidth', 2, 'color', col);
    end
    v0 = ctransform(c, x,u,y, 0);
    plot(y,v0, 'b.', 'MarkerSize', 25);
    set(gca, 'XTick', [], 'YTick', []);
    box on; axis tight; axis([0 1 min(v0)-m max(v0)+m]);
    SetAR(1/2);
    saveas(gcf, [rep 'c-transf-v-' num2str(round(10*r)) '.eps'], 'epsc');

    clf; hold on;
    for ieps=1:length(eps_list)
        epsilon = eps_list(ieps);
        mc = (ieps-1)/(length(eps_list)-1);
        col = [1-mc mc 0];
        %
        v = ctransform(c, x,u,y, epsilon);
        vg = ctransform(c, x,u,gy, epsilon);
        %
        uu = ctransform(c, y,v,x, epsilon);
        uug = ctransform(c, y,v,gx, epsilon);
        %
        vv = ctransform(c, x,uu,y, epsilon);
        vvg = ctransform(c, x,uu,gy, epsilon); % should match vg
        %
        plot(gx,uug, '-', 'LineWidth', 2, 'color', col);
    end 
    plot(x,u, 'r.', 'MarkerSize', 25);
    set(gca, 'XTick', [], 'YTick', []);
    box on; axis tight; axis([0 1 0-m 1+m]);
    SetAR(1/2);
    saveas(gcf, [rep 'c-transf-u-' num2str(round(10*r)) '.eps'], 'epsc');
end
