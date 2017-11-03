%%
% Test for W1 on graph using linprog.


addpath('../toolbox/');
rep = '../results/w1-graph/';
[~,~] = mkdir(rep);

% helpers
normalize = @(x)x/sum(x(:));

% generate a planar graph
N = 80;
X = randn(2,N);
X(2,:) = X(2,:)*.75;
T = delaunay(X(1,:),X(2,:));

I = [T(:,1);T(:,2);T(:,2);T(:,3);T(:,3);T(:,1)];
J = [T(:,2);T(:,1);T(:,3);T(:,2);T(:,1);T(:,3)];
[I,J] = deal(I(I<J), J(I<J));
P = length(I);

% gradient operator
grad = sparse( [1:P 1:P], [I;J]', [ones(P,1) -ones(P,1)], P,N);

clf;
triplot(T,X(1,:),X(2,:));


% edge length
w = sqrt( (X(1,I)-X(1,J)).^2 + (X(2,I)-X(2,J)).^2 );
w = w(:);

% source/sink
[~,s] = min(X(1,:));
d = sqrt( (X(1,s)-X(1,:)).^2 + (X(2,s)-X(2,:)).^2 );
[~,t] = max(d);

% target measures
test = 'diracs';
test = 'spread'; Q = 14; % Q-NN

switch test
    case 'diracs'
        mu = zeros(N,1); mu(s)=1;
        nu = zeros(N,1); nu(t)=1;
    case 'spread'
        v = [s t]; Mu = {};
        for i=1:2
            % distance to target
            x = X(:,v(i));
            d = sqrt( (x(1)-X(1,:)).^2 + (x(2)-X(2,:)).^2 );
            [d,Is] = sort(d, 'ascend'); d = d(1:Q); Is = Is(1:Q);
            Mu{i} = zeros(N,1);
            Mu{i}(Is) = normalize(sqrt( 1-d/d(end) ));
        end
        mu = Mu{1}; nu = Mu{2};
end

%%
% display the measures.

clf; hold on;
triplot(T,X(1,:),X(2,:), 'k');
% vertex according to mass
for i=1:N
    plot(X(1,i),X(2,i), '.', 'MarkerSize', .01 + 50*mu(i)/max(mu), 'color', 'r');
    plot(X(1,i),X(2,i), '.', 'MarkerSize', .01 + 50*nu(i)/max(mu), 'color', 'b');
end
axis equal; axis off;
saveas(gcf, [rep 'inputs.png'], 'png');

%%
% Solve Primal and Dual

% max <u,mu-nu> s.t.  |grad(u)| <= w
cvx_solver sdpt3 % SeDuMi %
cvx_begin quiet % sdp quiet
cvx_precision high;
variable u(N,1); %  real;
% norm( grad*u, Inf ) <= 1;
abs( grad*u ) <= w;
maximize( sum( u .* (mu-nu) ) );
cvx_end

% min |w.*f|_1 s.t.   div(f)=mu-nu
cvx_solver sdpt3 % SeDuMi %
cvx_begin quiet % sdp quiet
cvx_precision high;
variable f(P,1); %  real;
(grad'*f) - (mu-nu)  == 0;
minimize( sum( abs(f).*w ) );
% maximize( sum( abs(f(:)).*w(:) ) );
cvx_end

% plot flow
K = find(abs(f)>1e-3);
clf; hold on;
triplot(T,X(1,:),X(2,:), 'k-');
for k=K'
    plot(X(1,[I(k) J(k)]),X(2,[I(k) J(k)]), 'k', 'LineWidth', 4);
end
% vertex according to mass
for i=1:N
    plot(X(1,i),X(2,i), '.', 'MarkerSize', .01 + 50*mu(i)/max(mu), 'color', 'r');
    plot(X(1,i),X(2,i), '.', 'MarkerSize', .01 + 50*nu(i)/max(mu), 'color', 'b');
end
axis equal; axis off;
saveas(gcf, [rep 'flow.png'], 'png');
saveas(gcf, [rep 'flow.eps'], 'epsc');

% plot distance function
K = find(abs(f)>1e-3);
D = u-min(u(:)); D = D/max(D(:));
h = min(abs( (grad*u) ./ w ), 1);
clf; hold on;
% triplot(T,X(1,:),X(2,:), 'k');
% edge according to saturation index
for k=1:P
    col = h(k) * [0 1 0] + (1-h(k)) * [0 0 0];
    plot(X(1,[I(k) J(k)]),X(2,[I(k) J(k)]), 'color',  col,'LineWidth', 2);
end
% shortest path
%for k=K'
%    plot(X(1,[I(k) J(k)]),X(2,[I(k) J(k)]), 'k:', 'LineWidth', 3);
%end
% vertex according to distance u
for i=1:N
    plot(X(1,i),X(2,i), '.', 'MarkerSize', 30, 'color', [D(i) 0 1-D(i)]);
end
plot(X(1,s),X(2,s), 'r.', 'MarkerSize', 40);
plot(X(1,t),X(2,t), 'b.', 'MarkerSize', 40);
axis equal; axis off;
saveas(gcf, [rep 'potential.png'], 'png');
saveas(gcf, [rep 'potential.epsc'], 'epsc');
