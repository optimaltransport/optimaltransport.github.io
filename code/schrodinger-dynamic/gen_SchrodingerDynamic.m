%%
% test to display Brownian bridge interpolation.

addpath('../toolbox/');
rep = '../results/schrodinger-dynamic/';
[~,~] = mkdir(rep);

% number of points
N = 6;

randn('state', 66);
X = randn(2,N);
Y = randn(2,N)+.5;

P = 256; % number of step for brownian bridge simulation
K = 50; % numbef of bridges


a = 0+0*1i;
b = 1+1*1i;
sigma = .05;
Bc = brownian_bridge(P,K,a,b,sigma);

clf;
plot_colored(Bc);

%%
% Run Sinkhorn

epsilon = 5;

c = distmat(X,Y).^2;
mu = ones(N,1)/N; nu = mu;
options.niter = 5000;
[u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);

%%
% Display results using Brownian Bridges

sigma = .2*sqrt(epsilon);
if epsilon<.01
    sigma=0;
end

M = 100; % budget of path
G = round(gamma*M);
[I,J,NP] = find(G);
clf; hold on;
for s=1:length(I)
	a = X(1,I(s)) + 1i * X(2,I(s));
	b = Y(1,J(s)) + 1i * Y(2,J(s));
    K = NP(s);
    % generate a bunch of bridges
    Bc = brownian_bridge(P,K,a,b,sigma);
    % plot
    plot_colored(Bc);
end
plot(X(1,:), X(2,:), 'r.', 'MarkerSize', 30);
plot(Y(1,:), Y(2,:), 'b.', 'MarkerSize', 30);
axis off; axis equal;
saveas(gcf, [rep 'schrodinger-dynamic-eps' num2str(round(1000*epsilon)) '.png'], 'png');
