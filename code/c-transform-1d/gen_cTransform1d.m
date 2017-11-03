%%
% Test for c-transform and Laguerre cells in 1D.

addpath('../toolbox/');

rep = 'results/c-transform-1d/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);


% # input points
N = 7; 
% Grid size
q = 512;
rng('default'); rng(646);

M = rand(1,N);

% grid
t = linspace(0,1,q)';

% potentials
u = 0.1*randn(N,1);

% exponent/scale for the transform
if not(exist('p'))
    p = 2;
end
rho = 2;
if not(exist('epsilon'))
    epsilon = .01;
end

%%
% study influence of p

epsilon = 0;
p_list = linspace(.5,2,13);
lw = 2;
clf; hold on;
Col = [];
for ip=1:length(p_list)
    p = p_list(ip);
    s = (ip-1)/(length(p_list)-1);
    col = [0,s,1-s];  Col(end+1,:) = col;  
    % v(x,y) = min_i c((x,y),m_i) - u_i
    D = abs( repmat(M(:), [1 q]) - repmat(t(:)', [N 1])  );
    S = rho * D.^p - repmat(u,[1 q]);
    %
    [v,I] = min(S);
    if epsilon>0
        S = S - repmat(v, [size(S,1) 1]);
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) ) + v;
    end
    %
    plot(t, v, 'LineWidth', lw, 'Color', col);
end
%
plot(M, -u, 'r.', 'MarkerSize', 25);
axis tight; box on;
colormap(Col); caxis([p_list(1),p_list(end)]);
set(gca, 'FontSize', 20);
colorbar;
SetAR(1/2);
saveas(gcf, [rep 'c-transform-p.png'], 'png');
saveas(gcf, [rep 'c-transform-p.eps'], 'epsc');


%%
% study influence of epsilon

p = 1;
eps_list = linspace(0,.3,13);
lw = 2;
clf; hold on;
Col = [];
for ieps=1:length(eps_list)
    epsilon = eps_list(ieps);
    s = (ieps-1)/(length(eps_list)-1);
    col = [s,s,1-s];  Col(end+1,:) = col;  
    % v(x,y) = min_i c((x,y),m_i) - u_i
    D = abs( repmat(M(:), [1 q]) - repmat(t(:)', [N 1])  );
    S = rho * D.^p - repmat(u,[1 q]);
    %
    [v,I] = min(S);
    if epsilon>0
        S = S - repmat(v, [size(S,1) 1]);
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) ) + v;
    end
    %
    plot(t, v, 'LineWidth', lw, 'Color', col);
end
%
plot(M, -u, 'r.', 'MarkerSize', 25);
axis tight; box on;
colormap(Col); caxis([eps_list(1),eps_list(end)]);
set(gca, 'FontSize', 20);
colorbar;
SetAR(1/2);
saveas(gcf, [rep 'c-transform-eps.png'], 'png');
saveas(gcf, [rep 'c-transform-eps.eps'], 'epsc');