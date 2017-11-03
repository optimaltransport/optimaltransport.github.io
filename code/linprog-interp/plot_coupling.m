function plot_coupling(P,A,f)

col = {'r' 'b'};
ms = 25;

lw = 6;

P = P/sum(P(:));
N = size(P);
% marginals
p = {sum(P,2), sum(P,1)'};
        
[i,j,v] = find(P);

axis ij;
hold on;
if not(isempty(f))
    imagesc(f); axis image; axis off; colormap gray(256);
end
for q=1:length(i)
    plot([A{1}(1,i(q)) A{2}(1,j(q))],[A{1}(2,i(q)) A{2}(2,j(q))],'k-', 'LineWidth', lw*v(q)*N(1));
end
for k=1:2
    for i=1:N(k)
    	plot(A{k}(1,i),A{k}(2,i), '.', 'color', col{k}, 'MarkerSize', N(k)*ms*p{k}(i));
    end
end
    % drawnow;


end