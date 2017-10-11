%%
% Test for c-transform and Laguerre cells in 2D.

addpath('../toolbox/');

rep = 'results/c-transform-2d/';
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

% # input points
N = 5; 
% Grid size
q = 512;
rng('default'); rng(646);

% inputs
if not(exist('M'))
    M = rand(2,N);
if 1
    clf; axis([0 1 0 1]);
    M = [];
    for i=1:N
        [a,b] = ginput(1);
        M(:,end+1) = [a;b];
    end
end
end

% grid
t = linspace(0,1,q);
[Y,X] = meshgrid(t,t);
XY = [X(:)';Y(:)'];

% potentials
u = 0.1*randn(N,1);

% exponent/scale for the transform
if not(exist('p'))
    p = 2;
end
rho = 2;
if not(exist('epsilon'))
    epsilon = .01;
end


%%
% study influence of p

epsilon = 0;
p_list = [.5 1 1.5 2 2.5];
clf; hold on;
for ip=1:length(p_list)
    p = p_list(ip);
    
    % v(x,y) = min_i c((x,y),m_i) - u_i
    D = distmat(M,XY);
    S = rho * D.^p - repmat(u,[1 q*q]);
    [v,I] = min(S);
    if epsilon>0
        S = S - repmat(v, [size(S,1) 1]);
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) ) + v;
    end
    v = reshape(v,[q q]);
    I = reshape(I,[q q]);    
    
    clf; hold on;
    imagesc(t,t,v');
    contour(t,t,v', 20, 'k');
    % plot input with weight
    a = rescale(u, .5,2);
    for i=1:N
        plot(M(1,i), M(2,i), 'r.', 'MarkerSize', a(i)*30);
    end
    % plot contours of Laguerre cells
    for i=1:N
        A = double(I==i);
        contour(t,t,A',[1/2 1/2], 'k', 'LineWidth', 2);
    end
    colormap parula(256);
    caxis([min(v(:)) max(v(:))]);
    axis equal; axis off;
    saveas(gcf, [rep 'c-transform-p' num2str(round(10*p))  '.png'], 'png');
end


%%
% study influence of p


p = 1;
eps_list = [1e-3 .01 .1 .3];

clf; hold on;
for ieps=1:length(eps_list)
    epsilon = eps_list(ieps);
    
    % v(x,y) = min_i c((x,y),m_i) - u_i
    D = distmat(M,XY);
    S = rho * D.^p - repmat(u,[1 q*q]);
    [v,I] = min(S);
    if epsilon>0
        S = S - repmat(v, [size(S,1) 1]);
        v = - epsilon * log( sum(  exp( -S/epsilon ) ) ) + v;
    end
    v = reshape(v,[q q]);
    I = reshape(I,[q q]);   
    
    % logistic model
    R = exp( -S/epsilon );
    R = R ./ repmat( sum(R), [N 1] );
    R = permute(reshape(R,[N q q]), [2 3 1]);

    % render image of smoothed cells   
    Col = [[1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [0 1 1]; [1 0 1]; [.5 .5 .5]]';
    A = zeros(q,q,3);
    for s=1:N
        A = A + repmat(R(:,:,s), [1 1 3]) .* ...
                repmat( reshape(Col(:,s),[1 1 3]), [q q 1] );
    end
    %
    clf; hold on;
    imagesc(t,t,A);
    contour(t,t,v, 20, 'k');
    % plot input with weight
    a = rescale(u, .5,2);
    for i=1:N
        plot(M(2,i), M(1,i), 'r.', 'MarkerSize', a(i)*30);
    end
    % plot contours of Laguerre cells
    %for i=1:N
    %    A = double(I==i);
    %    contour(t,t,A,[1/2 1/2], 'k', 'LineWidth', 2);
    %end
    colormap parula(256);
    caxis([min(v(:)) max(v(:))]);
    axis equal; axis off; drawnow;
    saveas(gcf, [rep 'c-transform-eps' num2str(round(1e2*epsilon))  '.png'], 'png');
end


% '-eps' num2str(round(100*epsilon))

