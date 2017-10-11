d=4;




X=[0,0; 1,0; 0,1; 1,1]';
%X=randn(n,d);

prange= 1:0.05:4
maxeig=zeros(length(prange),1);
for ii=1:length(prange)
    p=prange(ii);    
    D=pairwiseDistance(X).^(p/2);
    
    simuls=100;
    ee=zeros(simuls,1);
    uu=zeros(simuls,1);
    
    
    h=simplex_grid_index_all(3,4,35)/4;
    
    %h=h(:,1:18);
    
    n= size(h,2);
    
    P=eye(n)-1/n*ones(n);
    
    r=zeros(n);
    for i=1:n
        for j=i+1:n
            r(i,j)=(mexEMD(h(:,i),h(:,j),D));
        end
    end
    r=r+r';
    
    r=r.^(2/p) % recover D.^2 Since D = r^1/p we get the formula.
    
    u=real(eig(P*r*P));
    maxeig(ii)=max(u);
end

figure('color','white');
plot(prange,maxeig,'-kx','linewidth',3)
xlabel('$p$ parameter for Wasserstein','interpreter','latex')
ylabel('Max. Eig. of Centered Distance Matrix','interpreter','latex')
%ylabel('$\max \text{eig}(\mathbf{J}\mathbf{D}_p^2\mathbf{J})$ ','interpreter','latex')
set(gca,'fontsize',16);

