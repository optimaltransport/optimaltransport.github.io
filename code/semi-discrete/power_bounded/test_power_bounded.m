% A file to demonstrate the power_bounded function
% Author: Muhammad F. Kasim (University of Oxford)

% test == 1 for drawing an ordinary voronoi diagram
% test == 2 for power diagram with random small weights
% test == 3 to see how long does it take for larger number of sites with random small weights
test = 2;

if test < 3
    N = 50; % number of sites
    x = rand(N, 1);
    y = rand(N, 1); % x and y coordinate of the sites
    if test == 1
        wts = zeros(N, 1); % weights of each coordinate (all-zeros weights produce voronoi diagram)
    elseif test == 2
        wts = rand(N, 1) * 0.001;
    end
    crs = [0,0; 0,1; 1,1; 1,0]; % the bounding box in clockwise order

    % get the power diagram
    [V,C] = power_bounded(x,y, wts, crs);

    % draw the resulted power diagram
    plot(x, y, 'r.'); % plot the sites' points
    hold on;
    for i = [1:N]
        xpoly = V(C{i},1);
        ypoly = V(C{i},2);
        xpoly = [xpoly; xpoly(1)];
        ypoly = [ypoly; ypoly(1)];
        plot(xpoly, ypoly, 'b-');
    end
    xlim([min(crs(:,1)), max(crs(:,1))]);
    ylim([min(crs(:,2)), max(crs(:,2))]);
    hold off;

elseif test == 3
    Ns = [100, 1000, 10000, 100000, 1000000];
    times = [];
    for N = Ns
        x = rand(N, 1);
        y = rand(N, 1);
        wts = rand(N, 1) * 0.0001;
        crs = [0,0; 0,1; 1,1; 1,0];
        
        % time the power diagram
        tic;
        power_bounded(x, y, wts, crs);
        times(end+1) = toc;
    end
    disp(times);
    
    loglog(Ns, times, 'bo', 'LineWidth', 2);
    hold on;
    loglog(Ns, times, 'b-', 'LineWidth', 2);
    set(gca, 'FontSize', 14)
    xlabel('Number of sites');
    ylabel('Time required (s)');
end
