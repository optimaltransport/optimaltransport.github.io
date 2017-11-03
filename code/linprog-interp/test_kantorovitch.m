%%
% Test Kantorovitch relaxation.

addpath('../toolbox/');
addpath('../toolbox/mexEMD/');
rep = 'results/kantorovitch/';
if not(exist(rep))
    mkdir(rep);
end
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);


do_round = 0;
p_exp = 1.05; % exponent for transport

N = 6;
N = Inf;

name = 'paris.jpg';
f = imread(name);

clf; imagesc(f); axis image; axis off; colormap gray(256);
hold on;

col = {'r' 'b'};
ms = 50;

if N==6
    load n6-positions
elseif N==Inf
    clf; imagesc(f); axis image; axis off; colormap gray(256);
    hold on;
    A = {[] []};
    for k=1:2
        s = 0;
        while s<N
            s = s+1;
            [x,y,b] = ginput(1);
            if b~=1
                break;
            end
            plot(x,y,'.', 'color', col{k}, 'MarkerSize', ms);
            A{k}(:,end+1) = [x;y];
        end
        N = size(A{k},2);
    end
else
    [n,p] = size(f);
    for k=1:2
        A{k} = [rand(1,N)*p;rand(1,N)*n];
    end
end

normalize = @(x)x/sum(x(:));
pN = 100; % number of masses
if N<=3
    p = {};
    for k=1:2
        r = rescale(rand(N,1),.1,1);
        r = [0;cumsum(r)]; r = round(r/r(end)*pN);
        r = diff(r);
        p{k} = r;
    end
else
    pN = 1;
    p = {};
    for k=1:2
        p{k} = normalize( rescale(rand(N,1),.1,1) );
    end
end

clf; imagesc(f); axis image; axis off; colormap gray(256);
hold on;
for k=1:2
    for i=1:N
    	plot(A{k}(1,i),A{k}(2,i), '.', 'color', col{k}, 'MarkerSize', N*ms*p{k}(i)/pN);
    end
    if N<=10
    for i=1:N
        text(A{k}(1,i)+25,A{k}(2,i),num2str(i), 'FontSize', 24, 'FontWeight', 'bold', 'color', col{k});
    end
    end
end
if N<=10
    saveas(gcf, [rep 'n' num2str(N) '-input-seeds.png'], 'png');
end

% distance matrix
c = distmat(A{1},A{2}).^p_exp;
c = c/max(c(:))*35;
if do_round
    c = round(c);
end

% best
figure(1);
[cost,P] = mexEMD(p{1},p{2},c);
clf; plot_coupling(P,A,f);
saveas(gcf, [rep 'n' num2str(N(1)) '-best.png'], 'png');

% display squares
if N(1)<=20
    figure(2);
    clf; hold on;
    for i=1:N
        for j=1:N
            x = i;
            y = j;
            e = P(i,j)/max(P(:));
            rectangle('Position',[y-e/2,x-e/2,e,e], 'FaceColor', [1/2,0,1/2], 'EdgeColor', 'k', 'LineWidth', 2);
        end
    end
    axis equal; axis ij;
    q = .5;
    axis([1-q N+q 1-q N+q]);
    box on;
    set(gca, 'XTick', 1:N, 'YTick', 1:N);
    set(gca, 'FontSize', 20);
    saveas(gcf, [rep 'n' num2str(N(1)) '-coupling.png'], 'png');
end

if N(1)<=20
    %% Display as histograms
    vmax = max( [p{1}(:);p{2}(:)] );
    for k=1:2
        figure(2+k);
        clf;
        bar(p{k}, col{k});
        axis tight; axis([.5 N+.5 0 vmax*1.2]);
        SetAR(1/3);
        set(gca, 'FontSize', 20);
        saveas(gcf, [rep 'n' num2str(N(1)) '-hist-' num2str(k) '.eps'], 'epsc');
    end
    %% Display as couplings
    figure(5);
    t = linspace(0,1,N(1));
    B = { [t*0;t] [t*0+.5;t] };
    clf; plot_coupling(P,B,[]);
    axis equal; axis off;
	saveas(gcf, [rep 'n' num2str(N(1)) '-best-lines.eps'], 'epsc');
end
