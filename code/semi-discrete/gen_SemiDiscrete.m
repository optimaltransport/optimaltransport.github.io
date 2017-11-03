%%
% Test for power diagrams

addpath('../toolbox/');
addpath('power_bounded/');
addpath('power_diagrams/');

rep = '../results/semi-discrete/';
[~,~] = mkdir(rep);

% number of sites
N = 200;
% x and y coordinate of the sites
xy = rand(N, 2);
sigma = .1;
xy = sigma*randn(N, 2)+1/2;

% weights
w = ones(N,1);

% the bounding box in clockwise order
bb = [0,0; 0,1; 1,1; 1,0];
% get the power diagram
[V,C] = power_bounded(xy(:,1),xy(:,2), w, bb);
% [PD, PDinf] = powerDiagramWrapper(xy, w);
% draw the resulted power diagram
plot_power(xy,V,C,bb);

q = .01;
axis([-q 1+q -q 1+q]);

% target histogram
Atgt = ones(N,1)/N;

% optimization run
tau = .0004 * N; % gradient step size
niter = 200;
iter_disp = round( [1:8 niter/16 niter/8 niter/4 niter/2 niter] );
iter_disp = 1:2:niter;
err = []; kdisp = 0;
for i=1:niter
    progressbar(i,niter);
    [V,C] = power_bounded(xy(:,1),xy(:,2), w, bb);
    A = area_power(xy,V,C,bb);
    if intersect(i,iter_disp)
        kdisp = kdisp+1;
        clf; plot_power(xy,V,C,bb); axis([-q 1+q -q 1+q]);
        saveas(gcf, [rep 'iteration-' num2str(kdisp), '.eps'], 'epsc');
        drawnow;
    end
    A = A/sum(A); % be sure to be normalized ...
    wnew = w - tau * (A(:)-Atgt(:));
    err(i) = norm(w-wnew);
    w = wnew;
end

% error decay
plot(err); axis tight;

% plot final OT
clf; plot_power(xy,V,C,bb, 2);
saveas(gcf, [rep 'matching.eps'], 'epsc');
axis([-q 1+q -q 1+q]);
