%%
% Comparison of various dual norms (aka integral probability metrics).

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
rep = ['results/dual-norms/'];
[~,~] = mkdir(rep);


% random test
normalize = @(x)x/sum(x(:));
k = @(x,y)-abs(x-y);
nx = 3; ny = 5;
x = randn(nx,1);
y = randn(ny,1);
p = normalize(rand(nx,1));
q = normalize(rand(ny,1));
d = rkhs_norm(k,x,p,x,p);
d = rkhs_norm(k,x,p,y,q);

name = '1dirac-1dirac'; 
name = '1dirac-2dirac';

% delta_0 <-> delta_t
xmax = 3;
tlist = linspace(-xmax,xmax,257)';

switch name
    case '1dirac-1dirac'
        p = 1; q = 1;
        x = @(t)0; 
        y = @(t)t;
    case '1dirac-2dirac'
        p = 1; q = [1;1]/2;
        x = @(t)0; 
        y = @(t)[-t/2;t/2];
    otherwise
        error('Unknown');
end

D = []; lgd = {};
% energy distance
k = @(x,y)-abs(x-y);
f = @(t)sqrt( rkhs_norm(k, x(t),p, y(t),q) );
D(:,end+1) = arrayfun(f,tlist);
lgd{end+1} = 'Energy';
% Gaussian RKHS
sigma = .3;
k = @(x,y)exp( -(x-y).^2 / (2*sigma^2) );
f = @(t)sqrt( rkhs_norm(k, x(t),p, y(t),q) );
D(:,end+1) = arrayfun(f,tlist);
lgd{end+1} = 'Gauss';
% W1 distance
D(:,end+1) = abs(tlist);
lgd{end+1} = 'W_1';
% Flat distance
D(:,end+1) = min(abs(tlist),1);
lgd{end+1} = 'Flat';

%% display
clf;
plot(tlist, D, 'LineWidth', 2); 
legend(lgd, 'Location', 'NorthWest');
axis tight;
SetAR(1/2);
set(gca, 'FontSize', 20);
saveas(gcf, [rep name '.eps'], 'epsc');
