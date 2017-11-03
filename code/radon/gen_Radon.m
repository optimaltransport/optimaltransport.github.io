%%
% Test for barycenters using the sliced transform via Radon transform

addpath('../toolbox/');
addpath('../toolbox/img/');
rep = '../results/radon/';
[~,~] = mkdir(rep);

normalize = @(x)x/sum(x(:));

% load the shapes
names = {'disk', 'two_disks'};
n0 = 300;
a = {};
for i=1:2
    f = load_image(names{i},n0);
    f = sum(f,3);
    f = rescale(f)>.5;
    if f(1)==1
        f=1-f;
    end
    % zero padding
    L = 100;
    g = zeros(n0+2*L);
    g(L+1:end-L,L+1:end-L)=f;
    %
    a{i} = normalize(g);
end



% angles for the Radon transform
P = 360;
theta = (0:P-1)'/P*180;

myradon = @(x)radon(x,theta);
myiradon = @(x)iradon(x,theta,'linear', 'Hann');

Ra = myradon(a{1});

clf
imageplot(myiradon(Ra));

Q = size(Ra,1);
t = linspace(0,1,Q)';
myInv = @(x)interp1(x, t, t, 'spline');
myInv = @(x)(inverseFunc(x,1)-1)/(Q-1);

C = {}; R = {}; iC = {};
for i=1:2
    % radon transform
    R{i} = myradon(a{i});
    % cumulant
    vmin = 0; % 1e-6*max(R{i}(:)); % to regularize
    C{i} = cumsum(R{i}+vmin);
    C{i} = C{i} ./ repmat(C{i}(end,:), [Q 1]);
    % inverse cumulative
    for k=1:P
        iC{i}(:,k) = myInv(C{i}(:,k));
    end
end

% barycenters
S = 9;
Ct = []; Rt = [];
for s=1:S
    r = (s-1)/(S-1);
    iCt = (1-r)*iC{1} + r*iC{2};
    for k=1:P
        Ct(:,k) = myInv( iCt(:,k) );
        Rt(:,k) = diff([0;Ct(:,k)]);
    end
    if s==1
        Rt = R{1};
    elseif s==S
        Rt = R{end};
    end
    at = myiradon(Rt);
    at = max(at,0);
    % concentrate
    v = sort(at(:)); v = v(round(.98*end));
    at = min(at,v);
    %
    col = [1-r;0;r];
    A = [];
    for l=1:3
        A(:,:,l) = at * col(l) + (max(at(:))-at);
    end
    imwrite(rescale(A(L+1:end-L,L+1:end-L,:)), [rep 'bary-' num2str(s) '.png'], 'png');
    A = [];
    for l=1:3
        A(:,:,l) = Rt * col(l) + (max(Rt(:))-Rt);
    end
    J = round(n0/4); % cropping
    imwrite(rescale(A(J+L+1:end-L-J,:,:)), [rep 'radon-' num2str(s) '.png'], 'png');
end
