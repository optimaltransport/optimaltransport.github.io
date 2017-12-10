%%
% Test for stochastic gradient descent on the semi-dual. 

%%
% You might want to check
%       http://nbviewer.jupyter.org/github/gpeyre/numerical-tours/blob/master/matlab/ml_4_sgd.ipynb
% for an introduction to SGD. 

addpath('../toolbox/');

rep = '../results/semi-discrete-sgd/';
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
Y = rand(d,m);
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
    Q = 200^2;
    tZ = linspace(0,1,sqrt(Q));
    [YZ,XZ] = meshgrid(tZ,tZ);
    Z = [XZ(:)'; YZ(:)'];
end

% Exponent p for Wasserstein_p
if not(exist('p'))
    p = 2;
end
% regularization strength. epsilon=0 <=> no regularization
if not(exist('epsilon'))
    epsilon = .01;
end

%%
% SGD

% solver: without or with Cesaro averaging (helps to stabilize)
solver = 'sga'; % NOT implemented YET
solver = 'sgd';

% number of point in each mini-batch. Useful specially if you use GPU
% parallel computations.
n = 1; 

% step size schedule.
tau0 = 1; % initial step size of SGD
k0 = 5; % control length of warmup-phase
switch solver  % use a slower tau decay for sga
    case 'sgd'
        tau = @(k)tau0/(1+k/k0);
    case 'sga'
        tau = @(k)tau0/(1+sqrt(k/k0));
end

% number of iterations
niter = 1000;
niter_check = 1; % check energy from time to time

% random colors for display in 2D
Col = rand(3,m);

% initial discrete potential
v = zeros(m,1);
E = []; % keep track of the evolution of the energy
for k=1:niter
    % do not do this at home, we only do this for book-keeping the energy
    if niter_check==1 || mod(k,niter_check)==1
        [uZ,gZ] = ctransform(Y, v, Z, p, epsilon);
        E(end+1) = mean(uZ)-sum(v.*b);
        clf;
        if d==2
            subplot(2,1,1);
        end
        plot(1+(0:length(E)-1)*niter_check, E, 'LineWidth', 2); axis tight;
        if d==2
            A = zeros(sqrt(Q),sqrt(Q),3);
            for ic=1:3
                for j=1:m
                    A(:,:,ic) = A(:,:,ic)  + Col(ic,j) * reshape( gZ(j,:), [sqrt(Q),sqrt(Q)]);
                end
            end
            subplot(2,1,2);
            hold on;
            imagesc(tZ,tZ,permute(A, [2 1 3]));
            plot(Y(1,:),Y(2,:), 'r.', 'MarkerSize', 20);
            axis image; axis off; 
        end
        drawnow;
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








