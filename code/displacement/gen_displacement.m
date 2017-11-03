%%
% Display displacement interpolation of dense and discrete measures.

addpath('../toolbox/');
addpath('../toolbox/mexEMD/');
addpath('./toolbox/img/');

rep = ['results/displacement/'];
[~,~] = mkdir(rep);

% helpers
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
ms = 15;
myplot = @(x,y,col)plot(x,y, '.', 'MarkerSize', ms, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'LineWidth', 2);
myplotS = @(x,y,ms,col)plot(x,y, '.', 'MarkerSize', ms, 'MarkerEdgeColor', col, 'MarkerFaceColor', col, 'LineWidth', 1);


% load the shapes
names = {'disk', 'two_disks'};
n0 = 200;
xy = {};
for i=1:2
    f = load_image(names{i},n0);
    f = sum(f,3);
    f = rescale(f)>.5;
    if f(1)==1
        f=1-f;
    end
    I = find(f(:));
    [x,y] = ind2sub([n0 n0],I);
    xy{i} = [x(:) y(:)];
    xy{i} = (xy{i}-1)/(n0-1);
end
N = [size(xy{1},1), size(xy{2},1)];


%%
% Display the point clouds.

clf;
subplot(2,1,1);
myplot(xy{1}(:,1), xy{1}(:,2), 'r');
subplot(2,1,2);
myplot(xy{2}(:,1), xy{2}(:,2), 'b');

%%
% Solve the linprog of OT between the dense clouds.

a = ones(N(1),1) / N(1);
b = ones(N(2),1) / N(2);
C = distmat(xy{1}',xy{2}').^2;
[cost,gamma] = mexEMD(a,b,C);

%%
% Render using density estimator.

[I,J,gammaij] = find(gamma);
tlist = linspace(0,1,5);
for k=1:length(tlist)
    t=tlist(k);
    col = [1-t;0;t];
    Xt = (1-t)*xy{1}(I,:) + t*xy{2}(J,:);
    % render as image
    K = 8;
    n1 = n0*K;
    Xt1 = round( K*Xt*(n0-1) + 1);
    A = zeros(n1,n1);
    %        A(Xt1(:,1) + n1*(Xt1(:,2)-1)) = gammaij;
    for m=1:size(Xt1,1)
        A(Xt1(m,1) + n1*(Xt1(m,2)-1)) = A(Xt1(m,1) + n1*(Xt1(m,2)-1)) + ...
            gammaij(m);
    end
    nF = 128*2; % width of the convolution kernel
    s = K*2;
    g = exp(-(-nF:nF).^2 / (2*s^2) );
    B = conv2(conv2(A, g, 'same')', g, 'same')';
    %
    Br = [];
    for r=1:3
        Br(:,:,r) = B * col(r) + (1-B);
    end
    clf; imageplot(Br); drawnow;
    imwrite(rescale(Br), [rep 'interp-dens-' num2str(k) '.png']);
end


%%
% Render on a small cloud of discrete samples.

P = 300;
for i=1:2
    I = randperm(N(i)); I = I(1:P);
    xyS{i} = xy{i}(I,:);
end

CS = distmat(xyS{1}',xyS{2}').^2;
[cost,gammaS] = mexEMD(ones(P,1)/P,ones(P,1)/P,CS);

[I,J,gammaij] = find(gammaS);
tlist = linspace(0,1,5);
for k=1:length(tlist)
    t=tlist(k);
    Xt = (1-t)*xyS{1}(I,:) + t*xyS{2}(J,:);
    % plot as points
    clf;
    hold on;
    for i=1:length(gammaij)
        myplotS(Xt(i,1), Xt(i,2), 20, [1-t 0 t]);
    end
    axis([0 1 0 1]);
    % dummy points
    plot([0 0], [1 1], '.', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
    axis equal; axis square;
    axis off; axis ij;
    drawnow;
	saveas(gcf, [rep 'interp-points-' num2str(k) '.eps'], 'epsc');
end
