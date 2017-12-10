%%
% Test for evolution of sinkhorn solution with epsilon in 1D.

addpath('../toolbox/');
repeps = '../results/sinkhorn-evol/epsilon/';
[~,~] = mkdir(repeps);
repiter = '../results/sinkhorn-evol/iterations/';
[~,~] = mkdir(repiter);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

N = 256;
x = linspace(0,1,N)';
normalize = @(x)x/sum(x(:));
gauss = @(m,s)exp( -(x-m).^2/(2*s^2) );


vmin = .05;
mu = normalize(vmin + gauss(.3,.09) );
nu = normalize(vmin + .6*gauss(.6,.06) + .4*gauss(.8,.04) );

clf;
plot([mu nu], 'LineWidth', 2);


clf;
plot(nu/max(nu), 'r', 'LineWidth', 2); axis tight;
box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
axis([1 N 0 1.02]);
SetAR(1/4);
saveas(gcf, [repeps 'input-1.eps'], 'epsc');
clf;
plot(mu/max(mu), 'b', 'LineWidth', 2); axis tight;
box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
axis([1 N 0 1.02]);
SetAR(1/4);
saveas(gcf, [repeps 'input-2.eps'], 'epsc');


c = distmat(x',x').^2;
options.verb = 0;
options.tol = 1e-10;
options.niter = 1000*3;

%%
% Evolution with epsilon.

q = 50;
r = .5;
eps0 = 1;
eps1 = 1e-3;
r = exp( log(eps1/eps0)/(q-1) );
eps_list = eps0*r.^((0:q-1)');



for k=1:length(eps_list)
    epsilon = eps_list(k);
    [u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);
    
    %%
    % Render in color
    t = (k-1)/(length(eps_list)-1);
    c1 = [t 0 1-t]; 
    c1 = .9*[t 0 1-t]; 
    %
    c1 = (1-t)*[1 0 1] + t*[0 1 0];
    c2 = [1 1 1];
    A = zeros(N,N,3);
    U = gamma/max(gamma(:));
    for i=1:3
        A(:,:,i) = c1(i)*U + c2(i)*(1-U);
    end
    
    % colormap 
    s = linspace(0,1,256)';
    CM = s*c1 + (1-s)*c2;
    
    % save images
    clf;
    imageplot(A); drawnow;
    imwrite(A, [repeps 'evol-img-' num2str(k) '.png'], 'png');
    
    %%% level-sets display %%%
    clf; hold on;
    imagesc(gamma);
    contour(gamma, 20, 'k');
    caxis([min(gamma(:)) max(gamma(:))]);
    axis equal; axis off;
    colormap(CM);
    saveas(gcf, [repeps 'evol-levelsets-' num2str(k) '.png'], 'png');
    

    %%% copula display %%%
    CP = copula(gamma);
    clf; hold on;
    imagesc(CP);
    contour(CP, 20, 'k');
    caxis([min(CP(:)) max(CP(:))]);
    axis equal; axis off;
    colormap(CM);
    saveas(gcf, [repeps 'evol-copula-' num2str(k) '.png'], 'png');
    
    %%% 3D display %%%
    zs = 0;
    clf; hold on;
    surf(x,x,zs + U');
    shading interp;
    plot3(x*0,x,zs + nu/max(nu), 'r', 'LineWidth', 2);
    plot3(x,x*0+1,zs + mu/max(mu), 'b', 'LineWidth', 2);
    colormap(CM);
    axis tight; view(55,45);
    view(65,60); lighting phong; camlight left; 
    box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    saveas(gcf, [repeps 'evol-3d-' num2str(k) '.png'], 'png');
    
end

return;


%% 
% Evolution with #iter

epsilon = 1e-3;
iter_list = [0 1 5 10 20 100 1000];
iter_list = 0:50;

for k=1:length(iter_list)
    options.niter = iter_list(k);
    [u,v,gamma,Wprimal,Wdual,err]  = sinkhorn_log(mu,nu,c,epsilon,options);
    
    %%
    % Render in color
    t = (k-1)/(length(iter_list)-1);
    c1 = [t 0 1-t]; 
    c1 = .9*[t 0 1-t]; 
    
    c1 = (1-t)*[1 0 1] + t*[0 1 0];
    
    c2 = [1 1 1];
    A = zeros(N,N,3);
    U = gamma/max(gamma(:));
    for i=1:3
        A(:,:,i) = c1(i)*U + c2(i)*(1-U);
    end
    
    clf;
    imageplot(A); drawnow;
    imwrite(A, [repiter 'evol-img-' num2str(k) '.png'], 'png');
    
    s = linspace(0,1,256)';
    CM = s*c1 + (1-s)*c2;
    
    zs = 0;
    clf; hold on;
    surf(x,x,zs + U');
    shading interp;
    plot3(x*0,x,zs + nu/max(nu), 'r', 'LineWidth', 2);
    plot3(x,x*0+1,zs + mu/max(mu), 'b', 'LineWidth', 2);
    colormap(CM);
    axis tight; view(55,45);
    view(65,60); lighting phong; camlight left; 
    box on; set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    saveas(gcf, [repiter 'evol-3d-' num2str(k) '.png'], 'png');
end