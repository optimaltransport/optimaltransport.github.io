%%
% Test for the computation of W-barycenter of Gaussian.

addpath('../toolbox/');
addpath('toolbox_tensors/');


d = 2;
K = 4; % #bary
S = @(u)u*u';
normalize = @(x)x/sum(x(:));
% generation of tensor from angle/aniso/scale
tensor = @(t,r,s)tensor_mult(tensor_creation(r, t), tensor_diag(s,s) );


for k=1:K
    C{k} = S(randn(d));
end

aniso = .7;
aniso = .98;
C = { ...
	tensor(.2,aniso,1), ...
	tensor(1.8,aniso,1.1), ...
	tensor(-.8,aniso,.9), ...
	tensor(0,0,1) };


lambda = normalize(rand(K,1));

niter = 50;
[B,err] = gbary(C,lambda,niter);
clf;
plot(log10(err)); axis tight;

Q = 9;
t = linspace(0,1,Q);
[Y,X] = meshgrid(t,t); X = X(:); Y = Y(:);
Lambda = [ (1-X).*(1-Y), (1-X).*Y, X.*(1-Y), X.*Y ];

B = {};
for s=1:Q^2
    progressbar(s, Q^2);
    lambda = Lambda(s,:);
    [B{s},err] = gbary(C,lambda,niter);
end

Col = [ [1;0;0] [0;1;0] [0;0;1] [1;.6;.2] ];

% eigenstructure
a = []; b = []; theta = [];
for s=1:Q^2
    lambda = Lambda(s,:);
    [U,E] = eig(B{s}); E = diag(E);
    a(s) = E(1); b(s) = E(2);
    theta(s) = atan2(U(2,1), U(1,1));
end

sc = .8/Q * 1/max(a);
clf; hold on;
for s=1:Q^2
    lambda = Lambda(s,:);
    [U,E] = eig(B{s}); E = diag(E);
    theta = atan2(U(2,1), U(1,1));
    col = sum( Col .* repmat(lambda(:)', [3 1]), 2 );
    ellipse_fill(E(1)*sc,E(2)*sc,theta,X(s),Y(s),col(:)');
end
h = .8/Q;
axis([-h 1+h -h 1+h]);
axis off;

rep = 'results/';
[~,~] = mkdir(rep);
saveas(gcf, [rep 'tensor-interp-' num2str(round(100*aniso)) '.png'], 'png');
% saveas(gcf, [rep 'tensor-interp-' num2str(round(100*aniso)) '.eps'], 'epsc');
