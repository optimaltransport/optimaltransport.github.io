Nx=4;
Ny=6;

%r=randomPointinDSimplex(Nx)+.05;
r=[.3,.1,.2,.4]'
%r=ones(Nx,1)/Nx;
%r=r/sum(r);
%c=randomPointinDSimplex(Ny)+.05;
c=[.1,.1,.2,.2,.3,.1]'
%c=c/sum(c);

%c=ones(Ny,1)/Ny;
X=1*randn(2,Nx)+0; 

X(:,4)=[-1.5;-2]
X(:,3)=[-1;3.5]
X(:,2)=[0;2]
X(:,2)=[.5;3.1]


% Y=randn(2,Ny)+1.3;%1.5;
% Y(:,end)=[3;5.5];

Y= [1.2817    2.6623    2.9484    0.8507    1    5.50000;
    1.7608    1.7519   -0.7284    2.2    0.0240    4.000]
C=pairwiseDistance(X,Y);%.^.2;


[d,T]=mexEMD(r,c,C);

% Get dual solutions by solving the degenerate linear system
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

ab(1:Nx)=ab(1:Nx)-median(ab);
ab(Nx+1:end)=ab(Nx+1:end)+median(ab);


scaling_factor=2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color','white')

subaxis(1,2,1)
hold on
for i=1:length(I)
    A=[X(:,I(i)),Y(:,J(i))];
    line(A(1,:),A(2,:),'color',[.5,.5,.5],'linewidth',max(1,round(30*T(I(i),J(i)))));
end


scatter(X(1,:),X(2,:),sqrt(r)*scaling_factor,'ob','filled','markeredgecolor','b')
hold on

scatter(Y(1,:),Y(2,:),sqrt(c)*scaling_factor,'sr','filled','markeredgecolor','r')    
set(gca,'fontsize',16);
set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 
grid on;
box on;
%export_fig dual_plot_1.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(1,2,2)
colormap hot
scatter(X(1,:),X(2,:),sqrt(r)*scaling_factor,ab(1:Nx),'o','filled','markeredgecolor','b','LineWidth',2)
hold on
scatter(Y(1,:),Y(2,:),sqrt(c)*scaling_factor,ab(Nx+1:Nx+Ny),'s','filled','markeredgecolor','r','LineWidth',2)    
cc = colorbar;

set(cc,'Location','manual')
set(gca,'fontsize',16);
set(gca,'yticklabel',{[]}) 
set(gca,'xticklabel',{[]}) 
grid on;
box on;
export_fig dual_plot.pdf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %figure('color','white')
% 
% hold on
% for i=1:length(I)
%     A=[X(:,I(i)),Y(:,J(i))];
%     line(A(1,:),A(2,:),'color',[.5,.5,.5],'linewidth',max(1,round(30*T(I(i),J(i)))));
% end
% 
% 
% scatter(X(1,:),X(2,:),sqrt(r)*scaling_factor,ab(1:Nx),'o','filled','markeredgecolor','b','LineWidth',1)
% hold on
% scatter(Y(1,:),Y(2,:),sqrt(c)*scaling_factor,ab(Nx+1:Nx+Ny),'s','filled','markeredgecolor','r','LineWidth',1)    
% colormap hot
% cc = colorbar;
% set(cc,'Location','manual')
% set(gca,'fontsize',16);
% set(gca,'yticklabel',{[]}) 
% set(gca,'xticklabel',{[]}) 
% grid on;
% box on;
% export_fig dual_plot_3.pdf
% 
% 
% 
