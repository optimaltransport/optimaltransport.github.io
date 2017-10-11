Nx=4;
Ny=6;

r=randomPointinDSimplex(Nx)+.05;
r=[.1,.3,.4,.2]'
%r=ones(Nx,1)/Nx;
r=r/sum(r);
c=randomPointinDSimplex(Ny)+.05;
c=[.1,.1,.2,.2,.3,.1]'
c=c/sum(c);

%c=ones(Ny,1)/Ny;
X=.8*randn(2,Nx)+0; 
X(:,end)=[-1.5;-2]

Y=.7*randn(2,Ny)+1.3;%1.5;
Y(:,end)=[3;5.5];
C=pairwiseDistance(X,Y);%.^.2;

[d,T]=mexEMD(r,c,C);

[I,J]=ind2sub([Nx Ny],find(T>0));

A=zeros(length(I)+1,Nx+Ny);
b=zeros(length(I)+1,1);
for i=1:length(I)
    A(i,I(i))=1;
    A(i,J(i)+Nx)=1;
    b(i)=C(I(i),J(i));    
end
A(end,:)=1; b(end)=0 % set sum = 0....

ab=A\b;

ab(1:Nx)=ab(1:Nx)-median(ab)
ab(Nx+1:end)=ab(Nx+1:end)+median(ab)



figure()

scaling_factor=2000;

subaxis(1,2,1)
hold on
for i=1:length(I)
    A=[X(:,I(i)),Y(:,J(i))];
    line(A(1,:),A(2,:),'color','black','linewidth',max(1,round(15*T(I(i),J(i)))));
end


scatter(X(1,:),X(2,:),sqrt(r)*scaling_factor,'b','filled','markeredgecolor','k')
hold on

scatter(Y(1,:),Y(2,:),sqrt(c)*scaling_factor,'r','filled','markeredgecolor','k')    
set(gca,'fontsize',16);
%axis off
grid on;
set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 

subaxis(1,2,2)
colormap hot
scatter(X(1,:),X(2,:),1/Nx*scaling_factor,ab(1:Nx),'^','filled','markeredgecolor','k')
hold on
scatter(Y(1,:),Y(2,:),1/Ny*scaling_factor,ab(Nx+1:Nx+Ny),'v','filled','markeredgecolor','k')    
cc = colorbar;
set(cc,'Location','manual')
set(gca,'fontsize',16);
set(gca,'yticklabel',{[]}) 
set(gca,'xticklabel',{[]}) 
%axis off
grid on;



