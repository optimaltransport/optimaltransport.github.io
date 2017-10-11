eps_list = [.01 .05 .1 .5 1 3 5 10 1000];
eps_list = [.01 .1 .3 .5 .7 1 100];

for ieps=1:length(eps_list)
    epsilon = eps_list(ieps);
    test_rates_regul;
    Deps{ieps} = D;
end

%%
% Display.

Alpha = .5; % transparency
mycolor = @(k,K)[ 1-(k-1)/(K-1), 0, (k-1)/(K-1) ];

for id=1:length(dList)
    d = dList(id);
    clf; hold on;
    lgd = {};
    for ieps=1:length(eps_list)
        epsilon = eps_list(ieps);
        lgd{ieps} = ['\epsilon=' num2str(epsilon)];
        plot(hx(nList), hy(mean(Deps{ieps}{id},2)), '-', 'Color', mycolor(ieps,length(eps_list)), 'LineWidth', 2);
    end    
    axis tight;
    legend(lgd);
    SetAR(1/2); box on;
    legend(lgd, 'Location', 'SouthWest');
    saveas(gcf, [rep 'influ-epsilon-d' num2str(d) '.eps'], 'epsc');
end
