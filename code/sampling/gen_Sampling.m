%%
% Test for discrete 2D sampling.

N = 512;

rep = '../results/sampling/';
if not(exist(rep))
    mkdir(rep);
end

name = 'ring';
name = 'peaks';

% helpers
rho = .8;
t = linspace(0,1,256);
cm = @(c0,c1)(t').^(rho)*c0(:)' + (1-(t').^(rho))*c1(:)';

gaussian = @(m,x,y) exp( -(x-y).^2 );

t = (0:N-1)'/N;
[Y,X] = meshgrid(t,t);
switch name
    case 'peaks'
        a = peaks(N);
        a = abs(a);
        col = [0 0 1];
    case 'ring'
        r = 1*sqrt((X-.5).^2+(Y-.5).^2);
        a = abs(sin(r))^2./(r+eps);
        a = cos(10*r).^2.* exp(-14*r.^2.2);
        a = min(a,.65);
%        surf(a); shading interp;
        col = [0 0 1];
end
a = a/sum(a(:));

% discrete samples
K = 1000;

for K=[100 500 1000 4000 10000]
    I = randsample(N^2,K,true,a(:));
    x = X(I); y = Y(I);

    clf; hold on;
    %imagesc(t,t,a);
    %colormap(cm(col,[1 1 1]));
    ms = 15;
    plot(x,y, '.', 'MarkerSize', ms, 'color', col);
    axis equal;
    axis([0 1 0 1]);  box on;
    set(gca, 'XTick', [], 'YTick', [], 'box', 'on');
    saveas(gcf, [rep name '-disc-' num2str(K) '.eps'], 'epsc');
end

col_cont = [1/2 0 1/2];
clf; hold on;
imagesc(t,t,a);
colormap(cm(col_cont,[1 1 1]));
q = 8;
contour(t,t,a, linspace(0,max(a(:)),q), 'color', col_cont);
caxis([0 max(a(:))]);
% axis off;
set(gca, 'XTick', [], 'YTick', [], 'box', 'on');
axis equal;
saveas(gcf, [rep name '-dens.png'], 'png');

%% pixelize
col_euler = [1 0 0];
b = a/max(a(:));
for n = [8 16 32 64]
P = N/n;
A = zeros(N,N,3);
for i=1:n
    for j=1:n
        seli = (i-1)*P+1:i*P;
        selj = (j-1)*P+1:j*P;
        u = mean(mean(b(seli,selj)));
        for s=1:3
            A(seli,selj,s) = col_euler(s)*u + (1-u);
        end
    end
end
A = A/max(A(:));
imwrite(A, [rep name '-pix-n' num2str(n) '.png'], 'png');
end
