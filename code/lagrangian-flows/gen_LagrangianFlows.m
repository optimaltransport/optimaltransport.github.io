%%
% Test for Folker-Planck equation discretization using Lagrangian gradient
% flows.

addpath('../toolbox/');

rep = '../results/lagrangian-flows/';
[~,~] = mkdir(rep);

% number of particles
n = 1000;
% number of neighbors for entropy estimation
k = 2;

% initial configuration
x = 2*rand(2,n)-1;

x = .4 * (2*rand(2,n)-1);
x(1,1:end/2) = x(1,1:end/2) + .5;
x(2,1:end/2) = x(2,1:end/2) + .5;
%
x(1,end/2+1:end) = x(1,end/2+1:end) -.3;
x(2,end/2+1:end) = x(2,end/2+1:end) - .5;


% potential function gradient
gP =  @(x)x;


% weight in front of entropy
kappa = 0.1*2;

mu = .01;

tau = .01;
niter = 150;
g = [];
r = 0;
for i=1:niter
    % NN-computations
    D = distmat(x,x);
    [Di,I] = sort(D, 1);
    I = I(2:k+1,:); Di = Di(2:k+1,:);
    % gradient for entropy  log(|x-xi|^2) -> (x-xi)/|x-xi|^2
    for s=1:2
        xs = x(s,:); xs = xs(:);
        h = repmat(x(s,:), [k 1]) - xs(I); h = h ./ (Di+mu).^2;
        g(s,:) = mean( h );
    end
    % gradient overall
    G = gP(x) - kappa * g;
    % gradient step
    x = x - tau*G;
    % display
    t = (i-1)/(niter-1);
    clf;
    plot(x(1,:), x(2,:), '.', 'Color', [t 0 1-t], 'MarkerSize', 20);
    axis([-1 1 -1 1]);
    axis equal; axis off;
    drawnow;
    if mod(i,20)==1
        r = r+1;
        saveas(gcf, [rep 'lagrange-flow-' num2str(r) '.eps'], 'epsc');
    end
end
