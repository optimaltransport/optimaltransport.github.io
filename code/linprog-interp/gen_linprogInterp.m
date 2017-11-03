%%
% test for linear programming and interpolation of discrete meaures.

addpath('../toolbox/');
addpath('../toolbox/mexEMD/');

test = 'weighted';
test = 'empirical';

rep = ['../results/linprog-interp/' test '/'];
[~,~] = mkdir(rep);

% helpers
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
myplot = @(x,y,ms,col)plot(x,y, 'o', 'MarkerSize', ms, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'LineWidth', 2);
myplot = @(x,y,ms,col)plot(x,y, 'o', 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 1);


% Dimensions  of the clouds.
n0 = 4000;
n1 = n0;

% Compute a first point cloud \(X_0\) that is Gaussian.
% and a second point cloud \(X_1\) that is Gaussian mixture.
randn('state', 666);
gauss = @(q,a,c)a*randn(2,q)+repmat(c(:), [1 q]);
X0 = randn(2,n0)*.3;
X1 = [gauss(n1/2,.5, [0 1.6]) gauss(n1/4,.3, [-1 -1]) gauss(n1/4,.3, [1 -1])];
% weights
normalize = @(a)a/sum(a(:));
switch test
    case 'weighted'
        p0 = normalize(rand(n0,1));
        p1 = normalize(rand(n1,1));
    case 'empirical'
        p0 = 2*ones(n0,1)/n0;
        p1 = 2*ones(n1,1)/n1;
end

%%
% Display the point clouds.
% The size of each dot is proportional to its probability density weight.

clf; hold on;
for i=1:length(p0)
    myplot(X0(1,i), X0(2,i), p0(i)*length(p0)*10, 'b');
end
for i=1:length(p1)
    myplot(X1(1,i), X1(2,i), p1(i)*length(p1)*10, 'r');
end
axis([min(X1(1,:)) max(X1(1,:)) min(X1(2,:)) max(X1(2,:))]); axis off;

%%
% Compute the cost matrix

C = repmat( sum(X0.^2)', [1 n1] ) + ...
    repmat( sum(X1.^2), [n0 1] ) - 2*X0'*X1;

%%
% Solve the linprog of OT

[cost,gamma] = mexEMD(p0,p1,C);

%%
% Compute displacement interpolation.

[I,J,gammaij] = find(gamma);
tlist = linspace(0,1,6);
clf;
for k=1:length(tlist)
    t=tlist(k);
    Xt = (1-t)*X0(:,I) + t*X1(:,J);
    % subplot(2,3,i);
    clf;
    hold on;
    for i=1:length(gammaij)
        myplot(Xt(1,i), Xt(2,i), gammaij(i)*length(gammaij)*6, [t 0 1-t]);
    end
    % title(['t=' num2str(t,2)]);
    axis([min(X1(1,:)) max(X1(1,:)) min(X1(2,:)) max(X1(2,:))]);
    % dummy points
    plot([min(X1(1,:)) max(X1(1,:))], [min(X1(2,:)) max(X1(2,:))], '.', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
    axis equal; axis square;
    axis off;
    drawnow;
    saveas(gcf, [rep 'interp-' num2str(k) '.eps'], 'epsc');
end
