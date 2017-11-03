%%
% Generate mixtures of Gaussian distributions


rot = @(t)[cos(t) sin(t);-sin(t) cos(t)];
cov = @(t,s)rot(-t)*diag(s)*rot(t);

% m=mean, C=cov
randncov = @(m,C,N)sqrtm(C)*randn(2,N)+repmat(m(:),[1 N]);
gaussian = @(m,C,X)1/sqrt(det(C))*exp( inv(C) );

C = cov(.5,[1 .1]);
m = [.5 1];
N = 1000;
X = randncov(m,C,N);

ms = 15;

clf;
plot(X(1,:), X(2,:), '.', 'MarkerSize', ms, 'color', [1 0 0]);
axis equal;