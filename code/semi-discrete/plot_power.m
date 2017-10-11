function plot_power(xy,V,C,bb, plot_mode)

N = size(xy,1);

ms = 25;
lw = 2;

if nargin<5
    plot_mode = 1;
end

% plot the sites' points
plot(xy(:,1), xy(:,2), 'r.', 'MarkerSize', ms);
hold on;
for i = [1:N]
    if not(isempty(C{i}))
        xpoly = V(C{i},1);
        ypoly = V(C{i},2);
        xpoly = [xpoly; xpoly(1)];
        ypoly = [ypoly; ypoly(1)];
        if plot_mode==1
            plot(xpoly, ypoly, 'b-', 'LineWidth', lw);
        else
            x1 = mean(xpoly); y1 = mean(ypoly);
            plot(x1, y1, 'b.', 'MarkerSize', ms);
            plot([xy(i,1) x1], [xy(i,2) y1], 'k-',  'LineWidth', lw);
        end
    end
end
xlim([min(bb(:,1)), max(bb(:,1))]);
ylim([min(bb(:,2)), max(bb(:,2))]);
hold off;
axis equal;
axis off;

end