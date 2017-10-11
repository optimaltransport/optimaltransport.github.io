function [cost] = plot_permutation(p,A,f)

col = {'r' 'b'};
ms = 40;

N = length(p);

axis ij;
hold on;
if not(isempty(f))
    imagesc(f); axis image; axis off; colormap gray(256);
end
for q=1:N
    plot([A{1}(1,q) A{2}(1,p(q))],[A{1}(2,q) A{2}(2,p(q))],'k-', 'LineWidth', 2);
end
for k=1:2
    plot(A{k}(1,:),A{k}(2,:),'.', 'color', col{k}, 'MarkerSize', ms);
end
% drawnow;


end