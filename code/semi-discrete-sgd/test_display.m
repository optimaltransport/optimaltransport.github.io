%%
% Test for stochastic gradient descent on the semi-dual -- Generates nice figures.  

%%
% You might want to check
%       http://nbviewer.jupyter.org/github/gpeyre/numerical-tours/blob/master/matlab/ml_4_sgd.ipynb
% for an introduction to SGD. 

addpath('../toolbox/');


% regularization strength. epsilon=0 <=> no regularization
if not(exist('epsilon'))
    epsilon = .01;
end


rep = ['../results/semi-discrete-sgd/eps-' num2str(round(1e3*epsilon)) '/'] ;
[~,~] = mkdir(rep);

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);

%%
% Parameters

% # points for the discrete measure nu
m = 10; 
% dimension of the space
d = 2;

%%
% Input measures

% function to draw n points from the continuous measure mu
X = @(n)rand(d,n);
% fixed points from the nu measure
if not(exist('Y'))
    % click and play
    Y = [];
    clf; hold on;
    while true
        axis([0 1 0 1]);
        [a,b,button] = ginput(1);
        plot(a,b, '.', 'MarkerSize', 20);
        if button==3
            break;
        end
        Y(:,end+1) = [a;b];
    end    
end
m = size(Y,2);
% weights on the nu measure
b = ones(m,1)/m;


%%
% WARNING: you should *not* do this in real life. We only do this for debug
%   purpose in small dimension, and geneate a dense sampling to monitor
%   evolution of the energy evaluated on a grid.

Q = 10000;
Z = X(Q);
if d==2
    % in dimension 2, use a uniform grid for simplicity
    Q = 400^2;
    tZ = linspace(0,1,sqrt(Q));
    [YZ,XZ] = meshgrid(tZ,tZ);
    Z = [XZ(:)'; YZ(:)'];
end

% Exponent p for Wasserstein_p
if not(exist('p'))
    p = 2;
end

%%
% SGD

% solver: without or with Cesaro averaging (helps to stabilize)
solver = 'sga'; % NOT implemented YET
solver = 'sgd';

% step size schedule.
tau0 = 1; % initial step size of SGD
k0 = 5; % control length of warmup-phase
switch solver  % use a slower tau decay for sga
    case 'sgd'
        tau = @(k)tau0/(1+k/k0);
    case 'sga'
        tau = @(k)tau0/(1+sqrt(k/k0));
end

% first run to get high quality value for the functional
niter = 10000;
n = 100; % number of point in each mini-batch.
v = zeros(m,1);
E = []; % keep track of the evolution of the energy
for k=1:niter
    progressbar(k,niter);
    % compute c-transform
    [u,g] = ctransform(Y, v, X(n), p, epsilon);
    % gradient of the semi-discrete energy, here
    G = b - mean(g,2);
    % SGD step
    v = v + tau(k) * G;
end
[uZ,gZ] = ctransform(Y, v, Z, p, epsilon);
E0 = mean(uZ)-sum(v.*b);

niter = 10000; % number of iterations
n = 1; % number of point in each mini-batch.
niter_check = 1; % check energy from time to time

% random colors for display in 2D
Col = distinguishable_colors(m)';

disp_list = [1:10:100 200:100:niter Inf];
kdisp = 1;

nrun = 4; % #runs
E = []; % keep track of the evolution of the energy
for irun=1:nrun    
    % initial discrete potential
    v = zeros(m,1);
    for k=1:niter
        progressbar(k,niter);
        % compute energy
        [uZ,gZ] = ctransform(Y, v, Z, p, epsilon);
        E(k,irun) = mean(uZ)-sum(v.*b);
        % do not do this at home, we only do this for book-keeping the energy
        if irun==1 && k==disp_list(kdisp)
            kdisp = kdisp+1;
            A = zeros(sqrt(Q),sqrt(Q),3);
            for ic=1:3
                for j=1:m
                    A(:,:,ic) = A(:,:,ic)  + Col(ic,j) * reshape( gZ(j,:), [sqrt(Q),sqrt(Q)]);
                end
            end
            clf; hold on;
            imagesc(tZ,tZ,permute(A, [2 1 3]));
            plot(Y(1,:),Y(2,:), 'r.', 'MarkerSize', 25);
            axis image; axis off;
            saveas(gcf, [rep 'iter-' num2str(k)], 'png');
        end
        % draw random samples
        x = X(n);
        % compute c-transform
        [u,g] = ctransform(Y, v, x, p, epsilon);
        % gradient of the semi-discrete energy, here
        %    E(v) = int v^c(x) da(x) - <v,b>
        %         ~ 1/n sum v^c(x_i) + <v,b>   [stochastic version]
        G = b - mean(g,2);
        % SGD step
        v = v + tau(k) * G;
    end
end



%% 
% display the energies

clf;
plot(1:niter, log(max(1e-5,1-E/E0)), 'LineWidth', 1); axis tight;
SetAR(1/2);
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'convergence.eps'], 'epsc');

