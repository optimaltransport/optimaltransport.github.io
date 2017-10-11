%%
% Visual display of convolutive kernels.

addpath('../toolbox/');
addpath('../toolbox/img/');
addpath('../toolbox/export_fig/');

rep = 'results/rkhs/';
[~,~] = mkdir(rep);


normalize = @(x)x/sum(x(:));

if not(exist('kernel'))
kernel = 'gaussian-medium';
kernel = 'gaussian-large';
kernel = 'gaussian-small';
kernel = 'energy-dist';
end

% grid points for display
q = 512;
% with padding
q1 = q*8;

% input 

rand('state', 1243);
N = 20;
mu = zeros(q,q);
I = randperm(q*q);
mu(I(1:N)) = 1/N;

nu = zeros(q,q);
I = randperm(q*q);
nu(I(1:N)) = 1/N;

mu = load_image('shape-1',q);
nu = load_image('shape-2',q);

mu = normalize(rescale(-mu));
nu = normalize(rescale(-nu));


xi = mu-nu;

% compute sqrt kernel 
u = 2 * [0:q1/2, -q1/2+1:-1]' / q1;
[Y,X] = meshgrid(u,u);
D = sqrt(X.^2+Y.^2);

% soft max
s=1/5; 
sm = @(t)s*log( 1+exp(t/s) );

switch kernel
    case 'energy-dist'
        r = .3;
        K = max(r-D,0);
        K = sm(r-D);
    case 'gaussian-small'
        sigma = .005;
        K = exp( -D.^2./(2*sigma^2) );
    case 'gaussian-medium'
        sigma = .02;
        K = exp( -D.^2./(2*sigma^2) );
    case 'gaussian-large'
        sigma = .05;
        K = exp( -D.^2./(2*sigma^2) );
end

Kh = fft2(K);
% Kh = max(real(Kh),0);

% for display
K1 = real( ifft2(sqrt(Kh)) );
K1 = fftshift(K1);
K1 = K1(end/2-q/2:end/2+q/2-1, end/2-q/2:end/2+q/2-1);

% zero padding
U = zeros(q1);
U(1:q,1:q) = xi;
xiC = real( ifft2( fft2(U) .* sqrt(Kh) ) );
xiC = xiC(1:q,1:q);

mu1 = mu/max(mu(:));
nu1 = nu/max(nu(:));
imagesc( cat(3, 2-nu1, 2-mu1-nu1, 2-mu1 )/2 );
axis image; axis off;
colormap jet(256);
% saveas(gcf, [rep 'input.png'], 'png');
export_fig([rep 'input.png'], '-m3');

disp_distrib(xiC);
axis image; axis off;
% saveas(gcf, [rep kernel '.png'], 'png');
export_fig([rep kernel '.png'], '-m3');

if strcmp(kernel, 'energy-dist')
    t = linspace(-1,1,q);
    [y,x] = meshgrid(t,t);
    r = .5*1e-3;
    A = -1 ./ sqrt(r + x.^2+y.^2); A = A/max(abs(A(:)));
    disp_distrib( A, 40 );
else
    disp_distrib(-K1, 40);
end
axis image; axis off;
% saveas(gcf, [rep kernel '-kernel.png'], 'png');
export_fig([rep kernel '-kernel.png'], '-m3');


