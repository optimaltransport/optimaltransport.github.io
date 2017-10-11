%%
% Plot of various divergence functionals

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
rep = ['results/divergences/'];
[~,~] = mkdir(rep);

t = linspace(0,6,1024)';

A = [ ...
    t .* log(t) - t + 1, ...
    abs(t-1), ...
    abs(sqrt(t)-1).^2, ...
    abs(t-1).^2, ...
];

lgd = {'KL', 'TV', 'Hellinger', '\chi^2'};

clf
plot(t, A, 'LineWidth', 2);
set(gca, 'FontSize', 20);
axis tight;
legend(lgd, 'Location', 'NorthWest');
SetAR(1/2);
saveas(gcf, [rep 'divergences.eps'], 'epsc');
