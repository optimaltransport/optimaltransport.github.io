%%
% Test for 1-D optimal transport.

N = 40; 
vmax = 20; 
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

rep = 'results/1d/';
[~,~] = mkdir(rep);

% generate Gaussian mixtures points in 1-D

x = 14 + 2*randn(N,1);
x = round(clamp(x,0,vmax));
y = 6 + 3*randn(N,1);
y = round(clamp(y,0,vmax));


ms = 30;
r = .2; E = 0*x + randn(N,1)*r;
% scatter plot
clf;
subplot(2,1,1);
plot(x, E, '.r', 'MarkerSize', ms); axis([0 vmax -.5 .5]); axis equal;
set(gca, 'YTick', []);
subplot(2,1,2);
plot(y, E, '.b', 'MarkerSize', ms); axis([0 vmax -.5 .5]); axis equal;
set(gca, 'YTick', []);


% histogram
subplot(2,1,1);
hx = hist(x,1:vmax); 
bar(1:vmax,hx, 'r'); axis tight;
subplot(2,1,2);
hy = hist(y,1:vmax); 
bar(1:vmax,hy, 'b'); axis tight;

%%
% animation on histograms

save_gif = 0;
K = 40; 
% sort
xs = sort(x);
ys = sort(y);
%
figure(2);
s = 0; F = [];
while true
for i=[1:K, K-1:-1:1]
    %
    t = (i-1)/(K-1);
    z = (1-t)*xs+t*ys;
    hz = hist(z,1:vmax); 
    clf;
    subplot(3,1,1);
    plot(z, E, '.', 'color', [1-t,0,t], 'MarkerSize', ms); axis([0 vmax -2 2]); % axis equal;
    set(gca, 'YTick', []);
    title('Optimal transport');
    subplot(3,1,2);
    bar(1:vmax,hz, 'FaceColor', [1-t,0,t]); axis([0,vmax,0,8]);
    title('Optimal transport histogram interpolation');
    subplot(3,1,3);
    bar(1:vmax,(1-t)*hx+t*hy, 'FaceColor', [1-t,0,t]); axis([0,vmax,0,8]);
    title('Linear histogram interpolation');
    drawnow;
    %
    if s<=K && save_gif
        s = s+1;
        f = getframe(gcf);
        F(:,:,:,s) = f.cdata;
    end
end
end


%%
% Save as .gif file.

if save_gif
    % find colormap
    A = permute(F, [1 2 4 3]);
    A = reshape(A, [size(A,1) size(A,2)*size(A,3) 3]);
    [~,map] = rgb2ind(uint8(A),254,'nodither');
    map(end+1,:) = 0; map(end+1,:) = 1;
    % convert
    im = [];
    for s=1:size(F,4);
        im(:,:,1,s) = rgb2ind(uint8(F(:,:,:,s)),map,'nodither');
    end
    % save
    imwrite(im+1,map,[rep 'interp-1d.gif'], ...
        'DelayTime',0,'LoopCount',inf);
end