%%
% Test for power diagrams

%%
% Test for power diagrams

addpath('power_bounded/');

rep = 'results/power-diagrams/';
[~,~] = mkdir(rep);

% number of sites
N = 50;
% x and y coordinate of the sites
sigma = .1;
xy = sigma*randn(N, 2)+1/2;
xy(1,:) = 1/2; % center points

% the bounding box in clockwise order
bb = [0,0; 0,1; 1,1; 1,0];

lambda = [1 5 10 100];
lambda = linspace(1,1.1,5);

for i=1:length(lambda)% weights
    w = ones(N,1);w(1) = lambda(i);    
    % get the power diagram
    [V,C] = power_bounded(xy(:,1),xy(:,2), w, bb);
    % draw the resulted power diagram
    plot_power(xy,V,C,bb);
    saveas(gcf, [rep 'power-' num2str(i), '.eps'], 'epsc');
    drawnow;
end
    