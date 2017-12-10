addpath('../toolbox/');
addpath('./toolbox/');
addpath('./shapes/'); 

%% Read data


N = 100; % size of images
epsilon = 1; % for N=100

N = 200; % size of images
epsilon = 2; % for N=200


names = {'annulus' 'cross' 'heart' '2disk'};
names = {'cat'  'star8'  'thinspiral'  'trefle'};

P = length(names);

H = zeros(N*N,P);
for i=1:P
    A = load_image(names{i}, N);
    A = sum(A,3);
    A = double(A>max(A(:))/2);
    if A(1)==1
        A = 1-A;
    end
    A = A + 1e-3;
    A = A / sum(A(:));
    H(:,i) = A(:);
end


rep = ['../results/barycenters-2d/' names{1}];
for i=2:P
    rep = [rep '-' names{i}];
end
rep = [rep '/'];
[~,~] = mkdir(rep);

n = size(H,1);
areaWeights = ones(n,1)/n;
H = H*n;

entropies = -sum(bsxfun(@times,H.*log(H),areaWeights),1);
maxEntropy = max(entropies);

%% Set up blur

imS = [N N];

if exist('imfilter') && 0
    % using image toolbox
    h = fspecial('gaussian',[1 max(imS)],epsilon);% hsize sigma
    h = h / sum(h);
    imBlur = @(x) imfilter(imfilter(x,h,'replicate'),h','replicate');
else
    blur = load_filtering('imgaussian', N);
    imBlur = @(x)blur(x,epsilon);
end

imagesc(imBlur(reshape(H(:,1),imS)));

blurColumn = @(x) reshape(imBlur(reshape(x,imS)),[],1);
blurAll2 = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));
blurAll = @(x) blurAll2(blurAll2(x));

% blurColumn = @(x) reshape(fastBlur(reshape(x,[targetSize,targetSize]),filterSize),[],1);
% blurAll = @(x) cell2mat(cellfun(blurColumn, num2cell(x,1), 'UniformOutput', false));

imagesc(reshape(blurAll(H(:,1)),imS));
axis equal;
axis off;

%% weights


% bilinear interpolation
q=5;
t = linspace(0,1,q);
[T,S] = meshgrid(t,t); S = S(:); T = T(:);
W = [(1-S).*(1-T) S.*(1-T) (1-S).*T S.*T]';
Q = size(W,2);

%% Compute barycenter

entropyLimit = Inf; % no limit

steps = linspace(0,1,5);% linspace(-.5,1.5,9);

averages = cell(length(steps),length(steps));
barycenters = cell(length(steps),length(steps));

options.niter = 200; %not 300
options.verb = 2;
options.tol = 1e-15;
resh = @(x)reshape(x, imS);
options.disp = @(x)imageplot(-resh(x));

col = [ [1 0 0]; [0 1 0]; [0 0 1]; [1 1 0]; [1/2 1 1/2]; [1/2 1/2 1] ];

B = {};
for i=1:q
    for j=1:q
        ij = (i-1)*q+j;
        w = W(:,ij)'; w = w/sum(w);
        %
        B{ij} = convolutionalBarycenter(H,w,areaWeights,blurAll,blurAll,entropyLimit,options);
        if isnan(sum(B{end}))
            warning('Problem');
            break;
        end
        u = rescale(resh(B{ij}));
        ui = u>1/2;
        % colorized version
        c = sum(col(1:P,:) .* repmat(w(:),[1 3]) );
        U = zeros(N,N,3); Ui = zeros(N,N,3);
        for r=1:3
            U(:,:,r) = u*c(r) + (1-u)*1;
            Ui(:,:,r) = ui*c(r) + (1-ui)*1;
        end
        clf; imageplot(U); drawnow;
        % save image
        imwrite(U, [rep 'measure-' num2str(i) '-' num2str(j) '.png'], 'png');
        imwrite(Ui, [rep 'shape-' num2str(i) '-' num2str(j) '.png'], 'png');
    end
end

