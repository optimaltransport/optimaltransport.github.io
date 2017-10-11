d=4;
n=35;




X=[0,0; 1,0; 0,1; 1,1]';
%X=randn(n,d);
p=2;
D=pairwiseDistance(X).^(p/2);

simuls=100;
ee=zeros(simuls,1);
uu=zeros(simuls,1);


jj=1;
%for jj = 1:simuls,
    
    h=simplex_grid_index_all(3,4,35)/4;
    %h=randomPointinDSimplex(d,n)';        
    
    
    ee(jj)=sum(-sum(h.*log(h)));
    disp(['Entropy : ',num2str(ee(jj))]);
    
    P=eye(n)-1/n*ones(n);
    
    r=zeros(n);
    for i=1:n
        for j=i:n
            r(i,j)=(mexEMD(h(:,i),h(:,j),D));%.^.1;
        end
    end
    r=r+r';
    
    u=real(eig(P*r*P));
    hist(u,1000)
    uu(jj)=max(u);
    disp(['Eigenval : ',num2str(uu(jj))]);
%end

scatter(ee,uu)