%%
% Test for sample complexity of sinkhorn loss.

addpath('../toolbox/');
rep = 'results/sample-complexity-regul/';
[~,~] = mkdir(rep);

% helps
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

options.niter = 400;
p = 1.4; % exponent
if not(exist('epsilon'))
    epsilon = 1;
end
sinkh_mode = 'mmd';


r = 15; % #n
q = 500; % #replications
%
r = 10; % #n
q = 400; % #replications

% use uniform distribution on a cube of dimension d
X = @(d,n)rand(d,n);

% log scale
nList = round( 10.^linspace(1,log10(300),r) );

dList = [2 5];
D = {};
for id=1:length(dList)
    d = dList(id);
    %    
    D{id} = [];
    for i=1:length(nList)
        progressbar(i,length(nList));
        n = nList(i);
        mu = ones(n,1)/n;
        for k=1:q
            D{id}(i,k) = sinkhorn_loss(X(d,n),X(d,n),mu,mu,epsilon,p, sinkh_mode, options);
        end
    end    
end


hx = @(x)log10(x);
hy = @(x)log10(x);
col = {'r' 'g' 'b' 'y'};

Alpha = .5; % transparency
clf; hold on; lgd = {};
for id=1:length(dList)
    plot(hx(nList), hy(mean(D{id},2)), '-', 'Color', col{id}, 'LineWidth', 2);
    lgd{id} = ['d=' num2str(dList(id))];
end
for id=1:length(dList)
%    plot(hx(nList), hy(D{id}), '.', 'Color', col{id}, 'MarkerSize', 10);
    shadedErrorBar(hx(nList), hy(mean(D{id},2)), std(hy(D{id}),[],2), {col{id} 'LineWidth' 2}, Alpha);
end
for id=1:length(dList)
    plot(hx(nList), hy(mean(D{id},2)), '-', 'Color', col{id}, 'LineWidth', 2);
end
axis tight;
legend(lgd);
SetAR(1/2); box on;
saveas(gca, [rep 'influ-d-eps' num2str(epsilon) '.eps'],'epsc');