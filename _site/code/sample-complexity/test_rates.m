%%
% Test for the convergence of empirical estimators.

addpath('../toolbox/');
addpath('../toolbox/toolbox-lsap/'); 

rep = 'results/sample-complexity/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

if not(exist('norm_type'))
norm_type = 'gaussian';
norm_type = 'energy';
norm_type = 'wass2';
end

r = 15; % #n
q = 500; % #replications

% linear scale
nList = round(linspace(10,3000,r));
% log scale
nList = round( 10.^linspace(1,3,r) );

dList = [2 3 5];
D = {};
for id=1:length(dList)
    d = dList(id);
    X = @(n)rand(d,n);
    s = 0*[1;zeros(d-1,1)]; % translation vector
    Y = @(n)X(n) + repmat(s,[1 n]);
    %    
    D{id} = [];
    for i=1:length(nList)
        progressbar(i,length(nList));
        n = nList(i);
        for k=1:q
            D{id}(i,k) = empirical_dist(norm_type, X(n), Y(n));
        end
    end    
end


Dtrue = 0;

hx = @(x)log10(x);
hy = @(x)log10(x);
col = {'r' 'g' 'b' 'y'};

Transp = .5;
clf; hold on; lgd = {};
for id=1:length(dList)
    plot(hx(nList), hy(mean(D{id},2)), '-', 'Color', col{id}, 'LineWidth', 2);
    lgd{id} = ['d=' num2str(dList(id))];
end
for id=1:length(dList)
%    plot(hx(nList), hy(D{id}), '.', 'Color', col{id}, 'MarkerSize', 10);
    shadedErrorBar(hx(nList), hy(mean(D{id},2)), 1*std(hy(D{id}),[],2), {col{id} 'LineWidth' 2}, Transp);
end
for id=1:length(dList)
    plot(hx(nList), hy(mean(D{id},2)), '-', 'Color', col{id}, 'LineWidth', 2);
end
axis tight;
legend(lgd);
SetAR(1/2); box on;






saveas(gcf, [rep norm_type '-rates.eps'], 'epsc');


return;

clf; % 
plot(nList, D, '.k', 'MarkerSize', 20); hold on;
loglog(nList, mean(D,2), '.-k', 'MarkerSize', 20, 'LineWidth', 2);
% axis([0 max(nList) 0 max(D(:))]);
axis tight;
SetAR(1/2); box on;


clf; hold on;
plot(nList, abs(D{1}-Dtrue), '.k', 'MarkerSize', 20);

S = sqrt( mean(D-Dtrue,2).^2 );
clf; hold on;
plot(nList, S, 'k', 'MarkerSize', 20);
